
# We use the Julia package "PowerModels" for implementation of the DC-OPF. 
# See https://lanl-ansi.github.io/PowerModels.jl/stable/ for more information
using PowerModels
PowerModels.silence()
using PGLib


# We use the test system "RTS-GMLC". The input file is in matpower format. 
# See https://github.com/GridMod/RTS-GMLC/blob/master/RTS_Data/FormattedData/MATPOWER/RTS_GMLC.m
# and https://matpower.org/docs/ref/matpower5.0/caseformat.html
# for more information
data = pglib((folder * "/data/RTS_GMLC.m"))

# SOME PRE-PROCESSING
for i in collect(keys(data["gen"]))
    # No unit commitment type constraints, i.e. no minimum power limits
    data["gen"][i]["pmin"] = data["gen"][i]["pmin"]*0

    # Remove the storage as time linking is not implemented
    if  data["gen"][i]["fuel"] == "Storage"
        data["gen"][i]["status"] = 0 
    end

end
# Ensure no angle difference limits
for branch in collect(keys(data["branch"]))
    data["branch"][branch]["angmin"] = -360
    data["branch"][branch]["angmax"] = 360
end


# Time series data from the RTS-GMLC
using DataFrames, CSV
include("functions.jl")


# Load time series data
RES_ts_data, load_ts_data, load_bus_data = load_time_series_data()
var_keys = variable_generator_keys(data, RES_ts_data)

# Load data for hydrogen production curve

# The curve is given for 1MW electrolyzer capacity and can be scaled linearly
df = CSV.read(folder * "/data/ELY_production_curve.csv", DataFrame)


