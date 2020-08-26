function parallel_oos!(id_year, nM, nQ, h, data, data_order, MNEMONIC, estim, ind_restr_states, nDraws, burnin, data_vintages,
                       data_vintages_year, unique_years, releases_per_year, oos_position, parfor_density_forecasts,
                       parfor_point_forecasts, parfor_rw_forecasts, parfor_outturn, parfor_parameters, parfor_states,
                       parfor_monthly_gdp, parfor_output_gap, parfor_BC_clean, parfor_EP_clean, parfor_BC, parfor_EP, parfor_T_INFL)

    @info("parallel_oos! > Starting $(unique_years[id_year]) out-of-sample");

    # --------------------------------------------
    # First release: estimate the coefficients
    # --------------------------------------------

    # Year ID
    v_year = 1;

    # Vintage ID (within each worker) <- for parfor_* output
    v = 1;

    # First vintage
    data_v = data_vintages_year[id_year][v_year];

    # Length data
    m = size(data)[1];

    # Last non-missing value per variable
    is_not_na   = ismissing.(data_v) .== false;
    last_not_na = zeros(1, nM+nQ);

    for var_id=1:nM+nQ
        last_not_na[var_id] = Int64(findall(is_not_na[:, var_id])[end]);
    end
    last_not_na = Array{Int64,2}(last_not_na);

    forecast_ends = last_not_na.+h;

    # Store: Random walk benchmark (must be before standardize_data! - or change rw_benchmark to use quarterly_position)
    # Output is in the right order
    parfor_rw_forecasts[:, :, v] = rw_benchmark(data_v, nM, nQ, data_order, last_not_na);

    # Generate outturn matrix (must be before standardize_data! - or change rw_benchmark to use quarterly_position)
    # Output must be reordered
    outturn_v = zeros(h, nM+nQ) |> Array{Union{Missing,Float64},2};
    for var_id=1:nM+nQ
        outturn_v[:, var_id] = data[last_not_na[var_id]+1:forecast_ends[var_id], var_id];
    end
    outturn_v[ismissing.(outturn_v)] .= NaN;

    # Standardize data_v
    data_v, _, quarterly_position, σʸ = standardize_data(data_v, nM, nQ, h, data_order, MNEMONIC);

    # SPF is unrestricted
    quarterly_position[MNEMONIC.=="GDP SPF"] .= 0.0;
    quarterly_position[MNEMONIC.=="INFL SPF"] .= 0.0;

    # Store: outturn
    parfor_outturn[:, :, v] = Array(outturn_v[:, data_order]);

    # Re-order last_not_na and forecast_ends
    last_not_na = last_not_na[data_order];
    forecast_ends = forecast_ends[data_order];

    # Run JuSSM
    distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, par, par_ind, par_size, distr_par =
    ssm_settings(data_v, h, nDraws, burnin, σʸ, quarterly_position, estim, ind_restr_states,
                 oos_details="$(unique_years[id_year]) out-of-sample");

    # Store: Model forecast
    for var_id=1:nM+nQ
        for draw=1:nDraws-burnin
            parfor_density_forecasts[:, var_id, draw, v] = distr_fcst[last_not_na[var_id]+1:forecast_ends[var_id], var_id, draw] .* σʸ[var_id];
        end
    end

    parfor_point_forecasts[:, :, v] = median(parfor_density_forecasts[:, :, :, v], dims=3);

    # Store: Parameters and states
    parfor_parameters[:,:]                  = chain_θ_bound;
    parfor_states[:, 1:size(distr_α)[2], :] = distr_α;


    # --------------------------------------------
    # Forecast
    # --------------------------------------------

    v += 1;

    no_vintages = length(data_vintages_year[id_year]);

    # This loop could be parallelised as well - for now we don't have enough computational power
    for v_year=2:no_vintages

        if (v_year==2) || (mod(v_year, 10)==0)
            @info("parallel_oos! ($(unique_years[id_year]) out-of-sample) > Vintage $v_year/$no_vintages");
        end

        data_v = data_vintages_year[id_year][v_year];

        # Last non-missing value per variable
        is_not_na   = ismissing.(data_v) .== false;
        last_not_na = zeros(1, nM+nQ);

        for var_id=1:nM+nQ
            last_not_na[var_id] = Int64(findall(is_not_na[:, var_id])[end]);
        end
        last_not_na = Array{Int64,2}(last_not_na);

        forecast_ends = last_not_na.+h;

        # Store: Random walk benchmark (must be before standardize_data! - or change rw_benchmark to use quarterly_position)
        # Output is in the right order
        parfor_rw_forecasts[:, :, v] = rw_benchmark(data_v, nM, nQ, data_order, last_not_na);

        # Generate outturn matrix (must be before standardize_data! - or change rw_benchmark to use quarterly_position)
        # Output must be reordered
        outturn_v = zeros(h, nM+nQ) |> Array{Union{Missing,Float64},2};
        for var_id=1:nM+nQ
            outturn_v[:, var_id] = data[last_not_na[var_id]+1:forecast_ends[var_id], var_id];
        end
        outturn_v[ismissing.(outturn_v)] .= NaN;

        # Standardize data_v
        data_v = data_v[:, data_order]./σʸ';

        # Store: outturn
        parfor_outturn[:, :, v] = Array(outturn_v[:, data_order]);

        # Re-order last_not_na and forecast_ends
        last_not_na = last_not_na[data_order];
        forecast_ends = forecast_ends[data_order];

        # Restrictions
        n_old = size(distr_par[1].y, 1)
        n_new = size(data_v, 2)

        if n_old > n_new
            data_v = [data_v missing .* ones(size(data_v, 1), n_old-n_new)];
            data_v[1:size(distr_par[1].y, 2), n_new+1:end] = distr_par[1].y[n_new+1:end, :]';
        end

        # Loop over variables and draws
        for draw=1:nDraws-burnin

            # Draw
            par_draw   = distr_par[draw];
            par_draw.y = data_v';
            α_draw, _  = kalman_diffuse!(par_draw, 0, 1, 1);

            # Store: Monthly GDP and output gap
            parfor_monthly_gdp[1:size(α_draw,2), draw, v] = (par_draw.Z[oos_position.GDP, 1:oos_position.last_GDP_state]' * α_draw[1:oos_position.last_GDP_state, :])' .* σʸ[oos_position.GDP];
            parfor_output_gap[1:size(α_draw,2), draw, v]  = (par_draw.Z[oos_position.GDP, 1:oos_position.last_GDP_state-1]' * α_draw[1:oos_position.last_GDP_state-1, :])' ./
                                                            (par_draw.Z[oos_position.GDP, oos_position.last_GDP_state] * α_draw[oos_position.last_GDP_state, :])';

            # Store: BC, EP and T_INFL
            parfor_BC_clean[1:size(α_draw,2), draw, v] = copy(α_draw[oos_position.BC,:]);
            parfor_EP_clean[1:size(α_draw,2), draw, v] = copy(α_draw[oos_position.EP,:]);
            parfor_BC[1:size(α_draw,2), draw, v]       = (par_draw.Z[oos_position.INFL, 1:oos_position.EP-1]' * α_draw[1:oos_position.EP-1, :])' .* σʸ[oos_position.INFL];
            parfor_EP[1:size(α_draw,2), draw, v]       = (par_draw.Z[oos_position.INFL, oos_position.EP:oos_position.T_INFL-1] * α_draw[oos_position.EP:oos_position.T_INFL-1, :])' .* σʸ[oos_position.INFL];
            parfor_T_INFL[1:size(α_draw,2), draw, v]   = (par_draw.Z[oos_position.INFL, oos_position.T_INFL] * α_draw[oos_position.T_INFL, :])' .* σʸ[oos_position.INFL];

            # Store: Model forecast
            for var_id=1:nM+nQ

                parfor_density_forecasts[:, var_id, draw, v] = (par_draw.Z[var_id, :]' * α_draw[:, last_not_na[var_id]+1:forecast_ends[var_id]])' .*σʸ[var_id];

            end
        end

        # Store point forecasts
        parfor_point_forecasts[:, :, v] = median(parfor_density_forecasts[:, :, :, v], dims=3);

        # Update v (within each worker)
        v += 1;
    end
end
