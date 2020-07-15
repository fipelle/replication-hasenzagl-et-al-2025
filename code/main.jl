#=
Name: main.jl
Description: Execution manager
Author: Filippo Pellegrino, f.pellegrino1@lse.ac.uk
=#

# Post burn-in size
n_distribution = nDraws-burnin;

# Load data
if run_type == 1
    data, date, nM, nQ, MNEMONIC = read_data(data_path);

# Load data for out-of-sample
else
    # Load general info
    MNEMONIC, nM, nQ, transf, transf_annex = read_data_info(data_info_path);

    # Load Fred data

    # - For the YoY% transformation
    if sum(transf .== 1) > 0
        df_fred_vintages = get_fred_vintages(fred_data_path, start_sample-Year(1), end_sample, oos_start_date);

    # - For the SPF transformation
    elseif sum(transf .== 1) == 0 && sum(transf .== 3) > 0
        df_fred_vintages = get_fred_vintages(fred_data_path, start_sample-Month(3), end_sample, oos_start_date);

    else
        df_fred_vintages = get_fred_vintages(fred_data_path, start_sample, end_sample, oos_start_date);
    end

    # Load local data
    df_local_vintages = get_local_vintages(local_data_path, end_sample, oos_start_date);

    # Merge vintages
    if size(df_fred_vintages,1) > 0 && size(df_local_vintages,1) > 0
        df_vintages = outerjoin(df_fred_vintages, df_local_vintages, on=[:date, :vintage_id]);

    # Use `df_fred_vintages` only
    elseif size(df_fred_vintages,1) > 0 && size(df_local_vintages,1) == 0
        df_vintages = copy(df_fred_vintages);

    # Use `df_local_vintages` only
    else
        df_vintages = copy(df_local_vintages);
    end

    # Array of vintages
    data_vintages, data_vintages_year, unique_years, releases_per_year = get_vintages(df_vintages, start_sample, MNEMONIC, transf, transf_annex, nM, h);

    # Last vintage
    data = data_vintages[end];
end

# Dimensions
m, n = size(data);


# ----------------------------------------------------------------------------------------------------------------------
# Execution: run_type == 1
# - Single iteration: it executes the code using the most updated datapoints
# ----------------------------------------------------------------------------------------------------------------------

if run_type == 1

    data, MNEMONIC, quarterly_position, σʸ = standardize_data(data, nM, nQ, h, data_order, MNEMONIC);
    data = [data; missing.*ones(h, nM+nQ)];

    # SPF is unrestricted
    quarterly_position[MNEMONIC.=="GDP SPF"] .= 0.0;
    quarterly_position[MNEMONIC.=="INFL SPF"] .= 0.0;

    @info("Data order: $MNEMONIC")

    # Run JuSSM
    distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, par, par_ind, par_size, distr_par =
        ssm_settings(data, h, nDraws, burnin, σʸ, quarterly_position, estim, ind_restr_states);

    # Transform variables to make them comparable with run_type==2
    data = data[1:end-h, :].*σʸ';

    # Save res in jld format
    save("./results/res$(res_name).jld", Dict("distr_α" => distr_α, "distr_fcst" => distr_fcst,
        "chain_θ_unb" => chain_θ_unb, "chain_θ_bound" => chain_θ_bound, "par" => par, "data_order" => data_order,
        "nDraws" => nDraws, "burnin" => burnin, "data" => data, "date" => date, "nM" => nM, "nQ" => nQ,
        "MNEMONIC" => MNEMONIC, "par_ind" => par_ind, "par_size" => par_size, "distr_par" => distr_par, "σʸ" => σʸ));


# ----------------------------------------------------------------------------------------------------------------------
# Execution: run_type == 2
# -  Out-of-sample
# ----------------------------------------------------------------------------------------------------------------------

elseif run_type == 2

    # --------------------------------------------
    # Initialise
    # --------------------------------------------

    # Quarterly position (this variable is temporary, as it is overwritten below)
    quarterly_position = [zeros(nM); ones(nQ)];

    if isempty(data_order) == false
        quarterly_position = quarterly_position[data_order];
        MNEMONIC = MNEMONIC[data_order];
    end

    # SPF is unrestricted
    quarterly_position[MNEMONIC.=="GDP SPF"] .= 0.0;
    quarterly_position[MNEMONIC.=="INFL SPF"] .= 0.0;

    @info("Data order: $MNEMONIC")

    # Dimensions
    k, size_θ       = ssm_settings(data, h, nDraws, burnin, ones(n), quarterly_position, estim, ind_restr_states, size_only_bool=true);
    n_vintages      = length(data_vintages);
    n_vintages_year = length(data_vintages_year);


    # --------------------------------------------
    # Run out-of-sample
    # --------------------------------------------

    @sync @distributed for id_year=1:length(unique_years)

        # Number of vintages per id_year
        n_vintages_id_year = length(data_vintages_year[id_year]);

        # Initialise output
        density_forecasts = Array{Float64}(undef, h, n, n_distribution, n_vintages_id_year);
        point_forecasts   = Array{Float64}(undef, h, n, n_vintages_id_year);
        rw_forecasts      = Array{Float64}(undef, h, n, n_vintages_id_year);
        outturn           = Array{Float64}(undef, h, n, n_vintages_id_year);
        parameters        = Array{Float64}(undef, size_θ, nDraws);               # full chain incl. burnin period
        states            = Array{Float64}(undef, k, m+h, n_distribution);

        # Monthly GDP and ouput gap
        monthly_gdp = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        output_gap  = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);

        # BC, EP and T_INFL
        BC_clean = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        EP_clean = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        BC       = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        EP       = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        T_INFL   = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);

        # Initialise to NaNs
        states      .= NaN;
        monthly_gdp .= NaN;
        output_gap  .= NaN;
        BC_clean    .= NaN;
        EP_clean    .= NaN;
        BC          .= NaN;
        EP          .= NaN;
        T_INFL      .= NaN;

        # Run parallel_oos!
        parallel_oos!(id_year, nM, nQ, h, data, data_order, MNEMONIC, estim, ind_restr_states, nDraws, burnin, data_vintages,
                      data_vintages_year, unique_years, releases_per_year, density_forecasts,
                      point_forecasts, rw_forecasts, outturn, parameters, states,
                      monthly_gdp, output_gap, BC_clean, EP_clean, BC, EP, T_INFL);

        # Save id_year-th chunk
        save("./results/res$(res_name)_chunk$(id_year).jld", Dict("id_year" => id_year,
             "density_forecasts" => density_forecasts,  "point_forecasts" => point_forecasts,
             "rw_forecasts" => rw_forecasts, "outturn" => outturn, "parameters" => parameters, "states" => states,
             "monthly_gdp" => monthly_gdp, "output_gap" => output_gap, "BC_clean" => BC_clean, "EP_clean" => EP_clean,
             "BC" => BC, "EP" => EP, "T_INFL" => T_INFL));
    end

    # Transform variables to make them comparable with run_type==1
    data     = data[:, data_order];
    MNEMONIC = MNEMONIC[data_order];

    # Save general settings as chunk0
    save("./results/res$(res_name)_chunk0.jld", Dict("estim" => estim, "data_vintages" => data_vintages,
       "data_vintages_year" => data_vintages_year, "unique_years" => unique_years,
       "releases_per_year" => releases_per_year, "nDraws" => nDraws, "burnin" => burnin, "data" => data, "date" => date,
       "nM" => nM, "nQ" => nQ, "MNEMONIC" => MNEMONIC, "data_order" => data_order));
end
