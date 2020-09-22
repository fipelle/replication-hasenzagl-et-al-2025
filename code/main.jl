#=
Name: main.jl
Description: Execution manager
=#

# Post burn-in size
n_distribution = nDraws-burnin;

# Load data

# ---------------------------------------------------------------------------------------------------------------------
# The following block of code is deprecated and it was used to load in-sample data from Excel files.
# ---------------------------------------------------------------------------------------------------------------------
#
# data, date, nM, nQ, MNEMONIC = read_data(data_path);
# ---------------------------------------------------------------------------------------------------------------------

# Load general info
MNEMONIC, nM, nQ, transf, transf_arg1, transf_arg2 = read_data_info(data_info_path);

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

# Re-order the data so that the quarterly series follow the monthly ones
df_vintages = df_vintages[!, vcat([:vintage_id, :date], Symbol.(MNEMONIC))];
sort!(df_vintages, :vintage_id);

# Array of vintages

# - For the YoY% transformation
if sum(transf .== 1) > 0
    data_vintages, data_vintages_year, unique_years, releases_per_year, unique_releases = get_vintages(df_vintages, start_sample-Year(1), end_sample, MNEMONIC, transf, transf_arg1, transf_arg2, nM, h);

# - For the SPF transformation
elseif sum(transf .== 1) == 0 && sum(transf .== 3) > 0
    data_vintages, data_vintages_year, unique_years, releases_per_year, unique_releases = get_vintages(df_vintages, start_sample-Month(3), end_sample, MNEMONIC, transf, transf_arg1, transf_arg2, nM, h);

else
    data_vintages, data_vintages_year, unique_years, releases_per_year, unique_releases = get_vintages(df_vintages, start_sample, end_sample, MNEMONIC, transf, transf_arg1, transf_arg2, nM, h);
end

# Last vintage
if run_type == 1 || run_type == 2

    # Find iis_release index
    diff_days = abs.(unique_releases .- iis_release);
    ind_iis_release = findall(diff_days .== minimum(diff_days))[1];

    # Select vintage and remove last h missings
    if run_type == 1
        @info("Selected vintage released on the $(unique_releases[ind_iis_release]) for in-sample estimation.");
    else
        @info("Selected vintage released on the $(unique_releases[ind_iis_release]) for conditional forecast estimation.");
    end

    data = data_vintages[ind_iis_release][1:end-h, :];

    # Generate corresponding vector of reference dates
    date = Dates.lastdayofmonth.(collect(range(start_sample,length=size(data,1),step=Dates.Month(1))));

elseif run_type == 3
    data = data_vintages[end];

else
    error("Wrong run_type!");
end

# Dimensions
m, n = size(data);

#=
------------------------------------------------------------------------------------------------------------------------
Run type == 1
------------------------------------------------------------------------------------------------------------------------
In-sample estimation: it executes the code using a single selected data vintage.
------------------------------------------------------------------------------------------------------------------------
=#
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

    # Remove the trailing h missing observations in data and the standardisation
    data = data[1:end-h, :].*σʸ';

    # Save res in jld format
    save("./results/res$(res_name).jld", Dict("distr_α" => distr_α, "distr_fcst" => distr_fcst,
        "chain_θ_unb" => chain_θ_unb, "chain_θ_bound" => chain_θ_bound, "par" => par, "data_order" => data_order,
        "nDraws" => nDraws, "burnin" => burnin, "data" => data, "date" => date, "nM" => nM, "nQ" => nQ,
        "MNEMONIC" => MNEMONIC, "par_ind" => par_ind, "par_size" => par_size, "distr_par" => distr_par, "σʸ" => σʸ));


#=
------------------------------------------------------------------------------------------------------------------------
Run type == 2
------------------------------------------------------------------------------------------------------------------------
Conditional forecast: it executes a series of conditional forecast on the basis of the in-sample coefficients,
and using a single selected data vintage. This option can be used only after having previously run the in-sample
estimation (run_type = 1).
------------------------------------------------------------------------------------------------------------------------
=#
elseif run_type == 2

    # Load parameters from in-sample output
    res_iis   = JLD.jldopen("$(pwd())/results/res$(res_iis_name).jld");
    σʸ        = permutedims(read(res_iis["σʸ"]));
    distr_par = read(res_iis["distr_par"]);
    nDraws    = read(res_iis["nDraws"]);
    burnin    = read(res_iis["burnin"]);
    @info("Loaded in-sample output.\\Replaced nDraws and burnin with values in $(res_iis_name).jld");

    # Standardize data with in-sample σʸ
    data, MNEMONIC, quarterly_position, _ = standardize_data(data, nM, nQ, h, data_order, MNEMONIC, σʸ);

    # Conditioning path
    if isempty(cond)
        data = [data; missing .* ones(h, n)];

    else
        cond_keys   = collect(keys(cond));
        cond_values = collect(values(cond));
        cond_pos    = vcat([findall(MNEMONIC .== j) for j in cond_keys]...);
        cond_h      = length.(cond_values);

        # Add h trailing missings if necessary
        data = [data; missing .* ones(maximum([h; cond_h]), n)];
        for j=1:length(cond_pos)
            data[m+1:m+cond_h[j], cond_pos[j]] = cond_values[j] ./ σʸ[cond_pos[j]];
        end
    end

    # SPF is unrestricted
    quarterly_position[MNEMONIC.=="GDP SPF"] .= 0.0;
    quarterly_position[MNEMONIC.=="INFL SPF"] .= 0.0;

    @info("Data order: $MNEMONIC")

    # --------------------------------------------
    # Compute conditional forecast output
    # --------------------------------------------

    # Initialise output
    # Note: at this stage the data is transposed wrt JuSSM_run.jl -> m, n = size(data) is correct.
    k               = size(distr_par[1].T)[1];
    m, n            = size(data);
    distr_cond_α    = zeros(k, m, nDraws-burnin);
    distr_cond_fcst = zeros(m, n, nDraws-burnin);

    # Loop over variables and draws
    for draw=1:nDraws-burnin

        # Draw
        par_draw   = distr_par[draw];
        par_draw.y = data';
        α_draw, _  = kalman_diffuse!(par_draw, 0, 1, 1);

        # Store states
        distr_cond_α[:, :, draw-burnin] = α_draw;

        # Store conditional forecasts (the predictions are for standardised data, as for the in-sample estimation)
        distr_cond_fcst[:, :, draw-burnin] = (par_draw.Z*α_draw)';
    end

    # Save res in jld format
    save("./results/res$(res_name).jld", Dict("distr_cond_α" => distr_cond_α, "distr_cond_fcst" => distr_cond_fcst));

#=
------------------------------------------------------------------------------------------------------------------------
Run type == 3
------------------------------------------------------------------------------------------------------------------------
Out-of-sample (real-time or pseudo, dependings on the settings in the Excel input).
------------------------------------------------------------------------------------------------------------------------
=#
elseif run_type == 3

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

        # Output gap and potential output
        output_gap  = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        potential_output = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);

        # BC, EP and T_INFL
        BC_clean = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        EP_clean = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        BC       = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        EP       = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);
        T_INFL   = Array{Float64}(undef, m+h, n_distribution, n_vintages_id_year);

        # Initialise to NaNs
        states      .= NaN;
        output_gap  .= NaN;
        potential_output .= NaN;
        BC_clean    .= NaN;
        EP_clean    .= NaN;
        BC          .= NaN;
        EP          .= NaN;
        T_INFL      .= NaN;

        # Run parallel_oos!
        parallel_oos!(id_year, nM, nQ, h, data, data_order, MNEMONIC, estim, ind_restr_states, nDraws, burnin, data_vintages,
                      data_vintages_year, unique_years, releases_per_year, oos_position, density_forecasts,
                      point_forecasts, rw_forecasts, outturn, parameters, states,
                      output_gap, potential_output, BC_clean, EP_clean, BC, EP, T_INFL);

        # Save id_year-th chunk
        save("./results/res$(res_name)_chunk$(id_year).jld", Dict("id_year" => id_year,
             "density_forecasts" => density_forecasts,  "point_forecasts" => point_forecasts,
             "rw_forecasts" => rw_forecasts, "outturn" => outturn, "parameters" => parameters, "states" => states,
             "output_gap" => output_gap, "potential_output" => potential_output,
             "BC_clean" => BC_clean, "EP_clean" => EP_clean, "BC" => BC, "EP" => EP, "T_INFL" => T_INFL));
    end

    # Adjust data order before saving chunk0
    data = data[:, data_order];

    # Save general settings as chunk0
    save("./results/res$(res_name)_chunk0.jld", Dict("estim" => estim, "data_vintages" => data_vintages,
       "data_vintages_year" => data_vintages_year, "unique_years" => unique_years,
       "releases_per_year" => releases_per_year, "nDraws" => nDraws, "burnin" => burnin, "data" => data,
       "df_vintages" => df_vintages, "start_sample" => start_sample, "end_sample" => end_sample,
        "transf" => transf, "transf_arg1" => transf_arg1, "transf_arg2" => transf_arg2,
       "nM" => nM, "nQ" => nQ, "MNEMONIC" => MNEMONIC, "data_order" => data_order));
end
