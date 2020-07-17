# ----------------------------------------------------------------------------------------------------------------------
# Initial settings
# ----------------------------------------------------------------------------------------------------------------------

# Load packages
using Distributed;
@everywhere using Dates;
@everywhere using FileIO;
@everywhere using LinearAlgebra;
@everywhere using Random;
@everywhere using Statistics;
@everywhere using FileIO, CSV, XLSX;
@everywhere using DataFrames;
@everywhere using Logging;
@everywhere using FredData;
@everywhere include("../../code/JuSSM/JuSSM.jl")
@everywhere using Main.JuSSM;

# General dependencies
@everywhere include("../../code/main_subroutines/data/read_data.jl");
@everywhere include("../../code/main_subroutines/data/data_vintages.jl");
@everywhere include("../../code/main_subroutines/data/standardize_data.jl");
@everywhere include("../../code/main_subroutines/oos/rw_benchmark.jl");
@everywhere include("../../code/main_subroutines/oos/parallel_oos!.jl");

# Model-specific dependencies
@everywhere include("../../code/ssm_settings_with_CBO.jl");

# Data paths
data_info_path = "./data/US_info.xlsx";
local_data_path = "./data/US_local.xlsx";
fred_data_path = "./data/US_fred.xlsx";

# Forecast horizon
h = 36; # forecast horizon [it is used when run_type is 1 or 2]


# ----------------------------------------------------------------------------------------------------------------------
# JuSSM settings
# ----------------------------------------------------------------------------------------------------------------------
# You can use two different settings for the initialisation and the execution
# ----------------------------------------------------------------------------------------------------------------------

nDraws = 10000;
burnin = 5000;


#=
------------------------------------------------------------------------------------------------------------------------
Run type
------------------------------------------------------------------------------------------------------------------------
1. Single iteration: it executes the code using the most updated datapoints
2. Out-of-sample (real-time or pseudo, dependings on the settings in the Excel input)
------------------------------------------------------------------------------------------------------------------------
=#

run_type = 2;
res_name = "";

# Out-of-sample options
oos_start_date = Dates.Date("01-01-2005", "dd-mm-yyyy");
start_sample = Dates.Date("01-01-1985", "dd-mm-yyyy");
end_sample = Dates.Date("30-06-2020", "dd-mm-yyyy");

#=
Data order is:
- CBO
- GDP
- GDP SPF
- URATE
- EMPL
- OIL
- INFL
- INFL SPF
- EXP INFL
=#
data_order = [6; 7; 8; collect(1:4); 9; 5];

# Estimate loading for the common trend (associated to EXP INFL)
estim = false;

# Linear restrictions in the filter and smoother (for the idiosyncratic trends when t=1)
ind_restr_states = []; # None


# ----------------------------------------------------------------------------------------------------------------------
# Execution
# ----------------------------------------------------------------------------------------------------------------------

# Random seed
Random.seed!(1);

# Run
include("../../code/main.jl");
