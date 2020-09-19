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
@everywhere using FileIO, CSV, JLD, XLSX;
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
h = 36;


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
1. In-sample estimation: it executes the code using a single selected data vintage.
2. Conditional forecast: it executes a series of conditional forecast on the basis of the in-sample coefficients,
   and using a single selected data vintage. This option can be used only after having previously run the in-sample
   estimation (run_type = 1).
3. Out-of-sample (real-time or pseudo, dependings on the settings in the Excel input).
------------------------------------------------------------------------------------------------------------------------
=#

run_type = 3;
res_name = "";
res_iis_name = ""; # used only in run_type == 2 to load the in-sample coefficients

#=
------------------------------------------------------------------------------------------------------------------------
Alfred settings
------------------------------------------------------------------------------------------------------------------------

1. iis_release and oos_start_date: "oos_start_date" is the date from which the code starts downloading the real-time
   vintages. In the in-sample estimation and when computing conditional forecasts, "oos_start_date" is used (for
   simplicity) to download a range of vintages from which only the data released closer to "iis_release" is selected
   and used. In other words, "iis_release" represents the date closer to the release of the vintage selected for
   in-sample estimation and to compute conditional forecasts. It is not used for the out-of-sample estimation.

2. start_sample and end_sample: first and last observations of interest.
------------------------------------------------------------------------------------------------------------------------
=#
iis_release = "";
oos_start_date = Dates.Date("01-01-2005", "dd-mm-yyyy");
start_sample = Dates.Date("01-01-1985", "dd-mm-yyyy");
end_sample = Dates.Date("31-08-2020", "dd-mm-yyyy");

# ----------------------------------------------------------------------------------------------------------------------
# Residual settings
# ----------------------------------------------------------------------------------------------------------------------

#=
Out-of-sample: position of the states and variables of interest [used only when run_type == 3]
- BC
- EP
- T_INFL
- GDP_trend
- GDP
- INFL
=#
oos_position = output_position(1, 5, 7, 10, 2, 7);

# Conditional forecast: conditioning path [used only when run_type == 2]
cond = Dict();

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
data_order = [7; 6; 8; collect(1:4); 9; 5];

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
