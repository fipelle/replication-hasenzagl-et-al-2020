# ----------------------------------------------------------------------------------------------------------------------
# Initial settings
# ----------------------------------------------------------------------------------------------------------------------

using Distributed;
include("code/read_data.jl");
include("code/tc_mwg.jl");
@everywhere include("code/Metropolis-Within-Gibbs/MetropolisWithinGibbs.jl")
@everywhere using DataFrames, Dates, FileIO, BSON, LinearAlgebra, Random, Statistics, XLSX;
@everywhere using Main.MetropolisWithinGibbs;

data_path = "./data/inflation.xlsx"; # Data file
end_presample_vec = [31, 12, 1998]; # End presample, day/month/year [it is used when run_type is 2 or 3]
h = 8; # forecast horizon [it is used when run_type is 1 or 3]


# ----------------------------------------------------------------------------------------------------------------------
# Metropolis-Within-Gibbs settings
# ----------------------------------------------------------------------------------------------------------------------

nDraws    = [60000; 40000]; # [number of draws in initialization; number of draws in execusion]
burnin    = nDraws .- 20000; # number of draws in the burn-in stage 
mwg_const = [0.025; 0.25]; # Initial constant. mwg_const might be adjusted to get an acceptance rate between 25% and 35%

#=
------------------------------------------------------------------------------------------------------------------------
Run type
------------------------------------------------------------------------------------------------------------------------
1. Single iteration: it executes the code using the most updated datapoints
2. Conditional forecast (you need to run option 1 first)
3. Out-of-sample: out-of-sample exercise, forecasting period starts after end_presample_vec
------------------------------------------------------------------------------------------------------------------------
=#

# In-sample run
run_type = 1;

if run_type == 1
	res_name = "iis"
elseif run_type == 2
	res_name = "cond"
elseif run_type == 3
	res_name = "oos"
end

#=
The variable `cond` is used only when run_type == 2.

It needs to be defined with a structure analogous to:

cond = [Dict("EMPL" => [51.30], "U" => [14.70], "OIL" => [16.61], "UOM" => [2.10]),
        Dict("EMPL" => [55.00], "U" => [10.00], "OIL" => [20.00]),
        Dict("U" => [15.00]),
        Dict("U" => [10.00]),
        Dict("U" => [15.00, 10.00])];
=#
cond = [];

# Load data
data, date, nM, nQ, MNEMONIC = read_data(data_path);


# ----------------------------------------------------------------------------------------------------------------------
# Execution
# ----------------------------------------------------------------------------------------------------------------------

# This random seed gives a chain similar to the one computed in Julia 0.6.2 for the paper
Random.seed!(2);

# Run code
include("code/tc_main.jl");
