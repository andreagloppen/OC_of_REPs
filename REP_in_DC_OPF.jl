
# Define current folder path - WITHOUT SLASH AT THE END

folder = @__DIR__
results_folder = folder * "/Results"

# Activate enviornment
using Pkg
Pkg.activate("env")

# Load input data and functions
include("input.jl")

# Activate remaining packages
using Gurobi, JSON3

# Define what time range to consider
T = range(1,2) #range(1,size(RES_ts_data)[1])

# REP specific inputs
allow_import = true # Is the REP allowed to import electricity? true or false
n_gen = 156         # RES to convert to an REP - idenficiation number
ely_cap = 1000       # Electrolyzer capacity at REP
lambda_h = 1.5      # Hydrogen price
steps = [1, 2000, 3000, 4000, 5000] # Electrolyzer curve steps - row numbers in the dataset for the hydrogen production curve (df)

# Network level modelling, "nodal", "zonal" or "copper" - nodal by default
network_level = "nodal"

for t in T 

    # Modify input data to current timestep
    # Changes time-varying input data 
    # This has to be done before adding the REP to get the REP cost function for the right RES available power
    data = mod_time_series_data(t, data)

    # Add the REP
    # Modifies the cost function of a renewable generation n_gen to represent the electrolyzer
    # Allows import by default, otherwise specify allow_import = false
    data = add_ely_general(n_gen, ely_cap, lambda_h, data, allow_import)


    # Modify network structure (network_level = ["nodal", "zonal", "copper"], "nodal" as standard)
    # If "nodal" there is no modification to the data. Otherwise, intra-regional or all tranmission line limits are relaxed to "M"
    data = change_transmission_cap(data, network_level)


    # Solve DC OPF
    # Using the PowerModels package implementation of a DC OPF
    # Select the solver - here using Gurobi, but Ipopt is also possible, see PowerModels documentation
    # Setting to output the LMPs of the system in the results dictionary
    result = solve_dc_opf(data, Gurobi.Optimizer, setting = Dict("output" => Dict("duals" => true)))



    # Additional result metrics
    result = calculate_gen_and_emissions(result, data) # Calculates total generation and emissions as well as emissions from each generator 
    result = calcucate_electrolyzer_setpoint(result, data, ely_cap, n_gen) # Calculates the power set-point of the electrolyzer in the REP
    result = calculate_curtailment(result, data, n_gen) # Calculates curtailment of the RES in the system - must be calculated after the electrolyzer set-point

    result = calculate_negative_cost_component(result, data, n_gen) # Calculates negative cost associated with electrolzyer consumption 
    

    # Save full result dictionary in JSON format - 
    open(results_folder * "/Dictionaries/results_t$(t).json", "w") do io JSON3.pretty(io, result, allow_inf=true)
    end
end



# Compute overall results in dataframes to save as csv files. 
status = DataFrame([col => [] for col in ["Status"; "Objective"; "Total_generation"; "Total_emission"; "Total_curtailment"; "Total_ely_consumption"; "Total_h2_production"; "Objective_negative"; "REP_cost"]])
lambdas = DataFrame([key => [] for key in collect(keys(data["bus"]))])
for t in T
    
    result = JSON3.read(results_folder * "/Dictionaries/results_t$(t).json", allow_inf=true, Union{Dict{String,Any},Array{Dict{String,Any}}})

    for key in collect(keys(result["solution"]["bus"]))
        push!(lambdas[!, key], result["solution"]["bus"][key]["lam_kcl_r"])
    end


    push!(status[!,"Status"], result["termination_status"]) 
    push!(status[!,"Objective"], result["objective"]) 
    push!(status[!,"Total_generation"], result["total_generation"]) 
    push!(status[!,"Total_emission"], result["total_emissions"]) 
    push!(status[!,"Total_curtailment"], result["total_curtailment"]) 
    push!(status[!,"Total_ely_consumption"], result["total_ely_consumption"]) 
    push!(status[!,"Total_h2_production"], result["total_h2_production"]) 
    push!(status[!,"Objective_negative"], result["objective_negative"]) 
    push!(status[!,"REP_cost"], result["solution"]["gen"]["$(n_gen)"]["pg_cost"]) 


end
CSV.write(results_folder * "/results_totals.csv", status)
CSV.write(results_folder * "/results_lambdas.csv", lambdas)

