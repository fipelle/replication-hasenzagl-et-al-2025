# ----------------------------------------------------------------------------------------------------------------------
# Initial settings
# ----------------------------------------------------------------------------------------------------------------------

using Distributed;

@everywhere using Dates;
@everywhere using FileIO;
@everywhere using LinearAlgebra;
@everywhere using Random;
@everywhere using Statistics;
@everywhere using FileIO, XLSX;
@everywhere using DataFrames;
@everywhere using Logging;

@everywhere include("../../code/main_subroutines/data/read_data.jl");
@everywhere include("../../code/main_subroutines/data/get_dataflow.jl");
@everywhere include("../../code/main_subroutines/data/standardize_data.jl");
@everywhere include("../../code/main_subroutines/oos/rw_benchmark.jl");
@everywhere include("../../code/main_subroutines/oos/parallel_oos!.jl");
@everywhere include("../../code/ssm_settings_no_core_CBO.jl");

@everywhere include("../../code/JuSSM/JuSSM.jl")
@everywhere using Main.JuSSM;

data_path = "./data/US.xlsx"; # Data file
h = 36; # forecast horizon [it is used when run_type is 1 or 3]


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
2. Out-of-sample: out-of-sample exercise
------------------------------------------------------------------------------------------------------------------------
=#

run_type = 1;

# when run_type = 2 this is the start date for the OOS
oos_start_date = Dates.DateTime("01-01-2005", "mm-dd-yyyy");

# when run_type == 2 it is an array of strings
res_name = "";

# used only when run_type == 2
cond = [];

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
