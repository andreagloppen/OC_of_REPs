

# Load time-series data
function load_time_series_data()
    # Imports all time seies data
    file_location = folder * "/data/timeseries_data_files/"

    # DA load in zone and distribution of zone loads onto nodes
    load_ts_data = CSV.read(file_location * "Load/DAY_AHEAD_regional_Load.csv", DataFrame)
    load_bus_data = CSV.read(file_location * "Load/BusData.csv", DataFrame)

    # All RES DA generation
    CSP = CSV.read(file_location * "CSP/DAY_AHEAD_Natural_Inflow.csv", DataFrame)
    hydro = CSV.read(file_location * "Hydro/DAY_AHEAD_hydro.csv", DataFrame)
    PV = CSV.read(file_location * "PV/DAY_AHEAD_pv.csv", DataFrame)
    RTPV = CSV.read(file_location * "RTPV/DAY_AHEAD_rtpv.csv", DataFrame)
    wind = CSV.read(file_location * "WIND/DAY_AHEAD_wind.csv", DataFrame)

    # Combinding all RES DA data into one CSV file
    RES_ts_data = innerjoin(CSP, hydro, PV, RTPV, wind, on = [:Year, :Month, :Day, :Period])

    # Return dataframes
    return RES_ts_data, load_ts_data, load_bus_data

end

# Get column key of RES generators
function variable_generator_keys(data, RES_ts_data)

    var_keys = []

    for i in collect(keys(data["gen"]))
        gen_name = data["gen"][i]["name"]
        try 
            test = RES_ts_data[!,gen_name]
            push!(var_keys,i)
        catch
        end
    end

    return var_keys

end


# Function to get data for a timestep t
function mod_time_series_data(t, data, var_keys = var_keys,  RES_ts_data = RES_ts_data, load_ts_data = load_ts_data, load_bus_data = load_bus_data)

    data_t = deepcopy(data)

    for i in var_keys
        gen_name = data["gen"][i]["name"]
        prod = RES_ts_data[t,gen_name]

        data_t["gen"][i]["pmax"] = prod/data["baseMVA"]

        # change maximum point on cost curve (just in case)
        pmax_cost = Int(data["gen"][i]["ncost"]*2-1)
        data_t["gen"][i]["cost"][pmax_cost] = data_t["gen"][i]["pmax"]

        if data_t["gen"][i]["pmax"] < data_t["gen"][i]["pmin"]
            println("Problem at key " * i)
            data_t["gen"][i]["gen_status"] = 0
        end
    end



    for i in collect(keys(data["load"]))

        # find the name of a load i
        load_bus = data["load"][i]["load_bus"]

        # find the row number of correspodning load name
        row_id = findfirst(bus -> bus == load_bus, load_bus_data[!,"Bus ID"])

        # find the share of total load for load i 
        proportion = load_bus_data[row_id,"Load Proportion"]

        # find the area code of load i
        area = load_bus_data[row_id,"Area"]

        # assign new load to load i, based on the share of area load and area load of time t 
        data_t["load"][i]["pd"] = proportion * load_ts_data[t,"$(area)"]/data["baseMVA"]
    
    end

    return data_t
end

# TODO check if more than a and b are needed from this function
function ely_approximated(df, steps)

    y = df.prod[steps]
    x = df.power[steps]

    I = range(1,length(x)-1)
    a = zeros(length(I))
    b = zeros(length(I))

    for i in I
        a[i] = (y[i] - y[i+1]) / (x[i] - x[i+1])
        b[i] = (y[i+1] - (a[i] * x[i+1]))
    end

    df.prod_approx = zeros(length(df.power))
    df.segm_approx = zeros(length(df.power))
    df.prod_m_approx = zeros(length(df.power))

    for j in range(1,length(df.power))
        for i in I
            if steps[i] <= j <= steps[i+1]
                df.prod_approx[j] = a[i] * df.power[j] + b[i]
                df.segm_approx[j] = i
            end

        end

        if j > 1
            df.prod_m_approx[j] = (df.prod_approx[j] - df.prod_approx[j-1])/(df.power[j] - df.power[j-1])
        end
    end

    # find unique values while ignoing the first element which is always zero
    prod_m_steps = unique(round.(df.prod_m_approx[2:5000], digits = 4)) 

    return a, b, prod_m_steps, df
end



function add_ely2(n_gen, steps, data, lambda_h, ely_cap, df = df)

    HPP_gen = deepcopy(data["gen"]["$(n_gen)"])

    # The electrolyzer capacity scaled after the possible full load hours
    C_E_FL = ely_capacity


    # Marginal production as approximated stepwise pieces
    a, b, prod_m_steps = ely_approximated(df, steps)
    power_steps = df.power[steps]

    # Assuming sell only 

    # Scale and flip for DA power
    p_da = C_E_FL .- power_steps * C_E_FL

    # The structure of the piecewise curve should be for increasing day-ahead power offer
    perm = sortperm(p_da)
    x = p_da[perm]

    a_sort = a[perm[2:length(perm)]]

    # We form the curve with increasing marginal cost, i.e. increasing slope.
    # We know the x-points and the slopes, we have the y points from h*lambda, find the b-s!

    point_1= zeros(length(x))
    point_2= zeros(length(x))
    b = zeros(length(x))

    for i in range(2,length(x)-1)
        point_1[i] = a_sort[i-1] * x[i] 
        point_2[i] = a_sort[i] * x[i] 
        b[i] = point_1[i] - point_2[i]
    end

    y = zeros(length(x))
    y[1] = a_sort[1] * x[1]
    for i in range(1,length(a))
        y[i+1] = a_sort[i] * x[i+1] + cumsum(b)[i]
    end



    # RES power at current time, scaled to MW
    p_RES = HPP_gen["pmax"] * data["baseMVA"]
    # Two different cases depending on if the RES power limits the electrolyzer consumption
    if p_RES <= C_E_FL

        # Should I limit the maximum point to the RES capacity?
        # I think this is ensured anyhow.

    elseif p_RES > C_E_FL

        # When there is more RES thwn ely cap, we sell some for free
        # We need two points to specify that the power is free for a range, not just a single point
        prepend!(x,[0,0])
        prepend!(y,[0,0])

        # How much is sold for free - the RES production that exceeds ely cap
        p_da_free = (p_RES - C_E_FL)

        # Push the x axis to the right to fit the zero offer piece
        x = x .+ p_da_free
        # Reset the first x point back to zero
        x[1] = 0


    end
    
    y = y .* lambda_h
    # The x-axis should be given in pu. Find it more intuitive to convert back at the very end here
    x = x ./data["baseMVA"]

    # structure the points into the same structure as the input file: mc_1 x_1 ... mc_n x_n
    temp = collect(Iterators.flatten(zip(x,y)))

    # make into an array
    HPP_cost = convert(Array{Float64,1}, temp)
    
    # change the cost parameters associated to the selected generator
    data["gen"]["$(n_gen)"]["cost"] = HPP_cost
    data["gen"]["$(n_gen)"]["ncost"] = length(x)
    data["gen"]["$(n_gen)"]["pmin"] = x[1]
    data["gen"]["$(n_gen)"]["pmax"] = x[end]
    data["gen"]["$(n_gen)"]["pres"] = p_RES/data["baseMVA"]

    # ensure that selected generators are online in input
    data["gen"]["$(n_gen)"]["gen_status"] = 1

    return data

end

function add_ely_import2(n_gen, steps, data, lambda_h, ely_cap, df = df)

    HPP_gen = deepcopy(data["gen"]["$(n_gen)"])

    # The electrolyzer capacity scaled after the possible full load hours
    #C_E_FL = ely_cap_dict[HPP_gen["name"]] 
    C_E_FL = ely_cap


    # Marginal production as approximated stepwise pieces
    a, b, prod_m_steps = ely_approximated(df, steps)

    prod_m_steps = a
    power_steps = df.power[steps]

    # Scale and flip for DA power
    p_da = C_E_FL .- power_steps * C_E_FL

    # The structure of the piecewise curve should be for increasing day-ahead power offer
    perm = sortperm(p_da)

    # 
    x = p_da[perm]

    #[0.0, 0.0, 0.951, 0.0]
    #[x_1, y_1,   x_2, y_2] where x is power and y is cost at the corresponding x

    # RES power at current time, scaled to MW
    p_RES = HPP_gen["pmax"] * data["baseMVA"]
    # Two different cases depending on if the RES power limits the electrolyzer consumption
    p_import_max = 0
    if p_RES <= C_E_FL

        # TODO I add the zero point but not sure if I actually need it. Compare with and without
        p_import_max = p_RES - C_E_FL

        x = x .+ p_import_max

        a_sort = a[perm[2:length(perm)]]

        # We form the curve with increasing marginal cost, i.e. increasing slope.
        # We know the x-points and the slopes, we have the y points from h*lambda, find the b-s!

        point_1= zeros(length(x))
        point_2= zeros(length(x))
        b = zeros(length(x))

        for i in range(2,length(x)-1)
            point_1[i] = a_sort[i-1] * x[i] 
            point_2[i] = a_sort[i] * x[i] 
            b[i] = point_1[i] - point_2[i]
        end

        y = zeros(length(x))
        y[1] = a_sort[1] * x[1]
        for i in range(1,length(a))
            y[i+1] = a_sort[i] * x[i+1] + cumsum(b)[i]
        end



    elseif p_RES > C_E_FL

        
        a_sort = a[perm[2:length(perm)]]

        # We form the curve with increasing marginal cost, i.e. increasing slope.
        # We know the x-points and the slopes, we have the y points from h*lambda, find the b-s!

        point_1= zeros(length(x))
        point_2= zeros(length(x))
        b = zeros(length(x))

        for i in range(2,length(x)-1)
            point_1[i] = a_sort[i-1] * x[i] 
            point_2[i] = a_sort[i] * x[i] 
            b[i] = point_1[i] - point_2[i]
        end

        y = zeros(length(x))
        y[1] = a_sort[1] * x[1]
        for i in range(1,length(a))
            y[i+1] = a_sort[i] * x[i+1] + cumsum(b)[i]
        end

        # When there is more RES thwn ely cap, we sell some for free
        # We need two points to specify that the power is free for a range, not just a single point
        prepend!(x,[0,0])
        prepend!(y,[0,0])

        # How much is sold for free - the RES production that exceeds ely cap
        p_da_free = (p_RES - C_E_FL)

        # Push the x axis to the right to fit the zero offer piece
        x = x .+ p_da_free
        # Reset the first x point back to zero
        x[1] = 0


    end
    
    y = y .* lambda_h
    # The x-axis should be given in pu. Find it more intuitive to convert back at the very end here
    x = x ./data["baseMVA"]

    # structure the points into the same structure as the input file: mc_1 x_1 ... mc_n x_n
    temp = collect(Iterators.flatten(zip(x,y)))

    # make into an array
    HPP_cost = convert(Array{Float64,1}, temp)
    
    # change the cost parameters associated to the selected generator
    data["gen"]["$(n_gen)"]["cost"] = HPP_cost
    data["gen"]["$(n_gen)"]["ncost"] = length(x)
    data["gen"]["$(n_gen)"]["pmin"] = x[1]
    data["gen"]["$(n_gen)"]["pmax"] = x[end]
    data["gen"]["$(n_gen)"]["pres"] = p_RES/data["baseMVA"]

    # ensure that selected generators are online in input
    data["gen"]["$(n_gen)"]["gen_status"] = 1

    return data

end


# Function to calculate emissions based on the result of the OPF
function calculate_gen_and_emissions(result, data)

    #store emissions values for each generator type
    gen_emissions_dict = Dict(
        "Oil"=>0.7434,
        "NG"=>0.6042,
        "Coal"=>0.9606,
        "Sync_Cond"=>0,
        "Nuclear"=>0,
        "Hydro"=>0,
        "Solar"=>0,
        "Wind"=>0,
        "Storage"=>0,
        "Ely"=>0

    )
    
    # Sum up the total emissions
    total_ems = 0
    total_gen = 0

    for gen in collect(keys(result["solution"]["gen"]))

        # Add a column to results for each generator
        result["solution"]["gen"][gen]["emissions"] = gen_emissions_dict[data["gen"][gen]["fuel"]]*result["solution"]["gen"][gen]["pg"]

        # Sum up over total
        total_ems = total_ems + result["solution"]["gen"][gen]["emissions"]
        total_gen = total_gen + result["solution"]["gen"][gen]["pg"]

    end

    # Add total emissions to result file
    result["total_generation"] = total_gen
    result["total_emissions"] = total_ems
    return result

end

# Function to calculate ely setpoint and h2 production based on the result of the OPF
function calcucate_electrolyzer_setpoint(result, data, ely_cap, n_gen, HPP_generators = [])

    if HPP_generators == []
        HPP_generators = ["$(n_gen)"]
    end
    # Sum up the total production/consumption
    total_ely = 0
    total_h = 0

    for gen in collect(keys(result["solution"]["gen"]))

        if gen in HPP_generators
            # Add a column of consumption for each electrolyzer
            result["solution"]["gen"][gen]["ely_set_point"] = data["gen"][gen]["pres"]-result["solution"]["gen"][gen]["pg"]

            # Adjust for the electrolyzer consumption
            if result["solution"]["gen"][gen]["ely_set_point"] > ely_cap/data["baseMVA"]
                result["solution"]["gen"][gen]["ely_set_point"]  = ely_cap/data["baseMVA"] #ely_cap_dict[data["gen"][gen]["name"]]
            end

            row_id = findfirst(x -> x >=  result["solution"]["gen"][gen]["ely_set_point"] ./ (ely_cap/data["baseMVA"]), df.power)

            result["solution"]["gen"][gen]["h2_production"] = df.prod[row_id] * ely_cap


        else
            result["solution"]["gen"][gen]["ely_set_point"] = 0
            result["solution"]["gen"][gen]["h2_production"] = 0
        end


        # Sum up total
        total_ely = total_ely + result["solution"]["gen"][gen]["ely_set_point"]
        total_h = total_h + result["solution"]["gen"][gen]["h2_production"]

    end

    # Add total column
    result["total_ely_consumption"] = total_ely
    result["total_h2_production"] = total_h
    return result

end


# DEPENDS ON ELY SETPOINT FUNCTION -  Function to calculate curtailment based on the result of the OPF
function calculate_curtailment(result, data, n_gen, HPP_generators = [])

    if HPP_generators == []
        HPP_generators = ["$(n_gen)"]
    end

    # Sum up the total curtailment
    total = 0
    for gen in var_keys

        # Add a column of curtailment for each generator
        result["solution"]["gen"][gen]["curtailment"] = data["gen"][gen]["pmax"]-result["solution"]["gen"][gen]["pg"]

        # Adjust for the electrolyzer consumption
        if gen in HPP_generators
            result["solution"]["gen"][gen]["curtailment"] = data["gen"][gen]["pres"] - result["solution"]["gen"][gen]["pg"]

            result["solution"]["gen"][gen]["curtailment"] = result["solution"]["gen"][gen]["curtailment"] - result["solution"]["gen"][gen]["ely_set_point"]
        end
        # Sum up total
        total = total + result["solution"]["gen"][gen]["curtailment"]
    end

    # Add total column
    result["total_curtailment"] = total
    return result

end

function calculate_negative_cost_component(result, data, n_gen)

    negative_component = 0
    for key in collect(keys(data["gen"]))
        result["solution"]["gen"][key]["pmax"] = data["gen"][key]["pmax"]
        if result["solution"]["gen"][key]["pg_cost"] < 0
            negative_component = negative_component + result["solution"]["gen"][key]["pg_cost"]
        end
    end
    result["objective_negative"] = negative_component

    return result
end

function additionality(HPP_generators = [], data = data, ely_cap_df = ely_cap_df)

    # Scale up RES source at a HPP location to have yearly matching with maximum ELY capacity
    for gen in HPP_generators
        row = filter(x -> x.node == parse(Int, gen), ely_cap_df).row_num
        res_capacity_factor = ely_cap_df.res_prod_year[row]/length(T_full)
        scale_factor = 1 + ely_cap_df.ely_cap[row][1]/res_capacity_factor[1]
        data["gen"][gen]["pmax"] = data["gen"][gen]["pmax"] * scale_factor
    end

    return data
end


function fixed_ely_consumption_2(data, node, HPP_gen, capacity, load_percentage)

    ely_no = parse(Int, HPP_gen)

    p_fixed = -capacity*load_percentage/data["baseMVA"]

    data["gen"][string(ely_no)] = Dict( 
    "ncost"      => 2,
    "model"      => 1,
    "gen_bus"    => node,
    "pmin"       => p_fixed,
    "pmax"       => 0,
    "mbase"      => 100.0,
    "source_id"  => Any["gen", ely_no],
    "index"      => ely_no,
    "cost"       => [p_fixed, -10000, 0, 0],
    "name"       => string(ely_no) * "_ELY",    
    "fuel"       => "Wind",
    "gen_status" => 1,
    "type"       => "ELY",
    "pg"         => 0.0,
    "apf"        => 0.0,
    "qc1max"     => 0.0,
    "shutdown"   => 0.0,
    "startup"    => 0.0,
    "qc2max"     => 0.0,
    "ramp_agc"   => 0.0,
    "qg"         => 0.0,
    "ramp_10"    => 0.0,
    "vg"         => 0.0,
    "pc2"        => 0.0,
    "qmax"       => 0.0,
    "qmin"       => 0.0,
    "qc1min"     => 0.0,
    "qc2min"     => 0.0,
    "pc1"        => 0.0,
    "ramp_q"     => 0.0,
    "ramp_30"    => 0.0
    )

    return data
end

function fixed_consumption_3(n_gen, fixed_load, data, ely_cap)
    p_RES = deepcopy(data["gen"]["$(n_gen)"]["pmax"]) # local production
    p_fixed = ely_cap*fixed_load/data["baseMVA"] # Electrolyzer fixed consumption

    A = 10000*data["baseMVA"] # MC that is higher than any other in the system

    if p_RES > p_fixed
        
        x = [0, p_RES-p_fixed, p_RES]
        y = [0, 0, A*(x[end]-x[end-1])]
    else
        x = [p_RES-p_fixed, 0]
        y = [A*(x[end-1]-x[end]),0]
    end

    temp = collect(Iterators.flatten(zip(x,y)))
    HPP_cost = convert(Array{Float64,1}, temp)

    # change the cost parameters associated to the selected generator
    data["gen"]["$(n_gen)"]["cost"] = HPP_cost
    data["gen"]["$(n_gen)"]["ncost"] = length(x)
    data["gen"]["$(n_gen)"]["pmin"] = x[1]
    data["gen"]["$(n_gen)"]["pmax"] = x[end]
    data["gen"]["$(n_gen)"]["pres"] = p_RES

    # ensure that selected generators are online in input
    data["gen"]["$(n_gen)"]["gen_status"] = 1
end


function add_ely_general(n_gen, ely_cap, lambda_h, data,  allow_import = true, hybrid = true, fixed_load = false, steps = steps, df = df)


    if hybrid == true


        if fixed_load == false

            if allow_import == false
                data = add_ely2(n_gen, steps, data, lambda_h, ely_cap)

            else
                data = add_ely_import2(n_gen, steps, data, lambda_h, ely_cap)
            
            end

        else
            data = fixed_consumption_3(n_gen, fixed_load, data, ely_cap)
        end



    else 

        gen_bus = data["gen"]["$(n_gen)"]["gen_bus"]

        new_gen = 159

        data["gen"]["$(new_gen)"] = deepcopy(data["gen"]["$(n_gen)"])
        data["gen"]["$(new_gen)"]["name"] = "$(gen_bus)_ELY_" * "$(new_gen)"
        data["gen"]["$(new_gen)"]["fuel"] = "Ely"
        data["gen"]["$(new_gen)"]["index"] = "$(new_gen)"
        data["gen"]["$(new_gen)"]["source_id"][2] = new_gen

        data["gen"]["$(new_gen)"]["pmax"] = 0.0

        
        data = add_ely_import2(new_gen, steps, data, lambda_h, ely_cap)

        if fixed_load != false
            p_fixed = -ely_cap*fixed_load/data["baseMVA"]
            data["gen"]["$(new_gen)"]["pmax"] = 0
            data["gen"]["$(new_gen)"]["pmin"] = p_fixed
            data["gen"]["$(new_gen)"]["ncost"]= 2
            data["gen"]["$(new_gen)"]["cost"] = [p_fixed, -100000, 0, 0]
        end


    end






    return data

end


function add_res(n_gen, scale_factor, data)
    data["gen"][n_gen]["pmax"] = data["gen"][n_gen]["pmax"] * scale_factor

end


function change_transmission_cap(data, network_level = "nodal", M = 400)

    for n_branch in collect(keys(data["branch"]))
        if network_level == "copper"
            data["branch"]["$(n_branch)"]["rate_a"] = M
            data["branch"]["$(n_branch)"]["rate_b"] = M
            data["branch"]["$(n_branch)"]["rate_c"] = M
        elseif network_level == "zonal"
            if last(digits(data["branch"]["$(n_branch)"]["t_bus"])) == last(digits(data["branch"]["$(n_branch)"]["f_bus"]))
                data["branch"]["$(n_branch)"]["rate_a"] = M
                data["branch"]["$(n_branch)"]["rate_b"] = M
                data["branch"]["$(n_branch)"]["rate_c"] = M

                data["branch"]["$(n_branch)"]["angmin"] = -360
                data["branch"]["$(n_branch)"]["angmax"] = 360

                data["branch"]["$(n_branch)"]["b_to"] = 1
                data["branch"]["$(n_branch)"]["b_fr"] = 1
                data["branch"]["$(n_branch)"]["br_r"] = 0.0001
                data["branch"]["$(n_branch)"]["br_x"] = 0.0001


            else

            end
        end
    end

    return data
end
