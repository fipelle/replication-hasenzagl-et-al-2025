# ----------------------------------------------------------------------------------------------------------------------
# Initial settings
# ----------------------------------------------------------------------------------------------------------------------

using JLD, FileIO;
using DataFrames;
using XLSX;
using Colors, PlotlyJS;
using ORCA, PlotlyBase;
using Dates, Statistics;

# Colors
c1 = "rgba(0, 48, 158, .75)";
c2 = "rgba(255, 0, 0, .75)";
c3 = "rgba(255, 190, 0, .75)";


# ----------------------------------------------------------------------------------------------------------------------
# Custom dependencies
# ----------------------------------------------------------------------------------------------------------------------

"""
    dropdims_median(X::Array{Float64,3})

Take median and remove unit dimension.
"""
dropdims_median(X::Array{Float64,3}) = dropdims(median(X, dims=2), dims=2);

"""
    load_oos_recon(nyears, model_folder, is_baseline; remove_forecast_path=true, vintages_to_exclude=0)

Load out-of-sample reconstruction results from the output chunks.
"""
function load_oos_recon(nyears, model_folder, is_baseline; remove_forecast_path=true, vintages_to_exclude=0)

    # Load chunk0
    chunk0 = jldopen("$(model_folder)/results/res_chunk0.jld");

    # Initial settings
    data_order        = read(chunk0["data_order"]);
    releases_per_year = read(chunk0["releases_per_year"]);

    # Dimensions
    h    = size(read(JLD.jldopen("$(model_folder)/results/res_chunk1.jld")["point_forecasts"]))[1];
    T, n = size(read(chunk0["data"]));

    # Number of releases to keep
    releases_to_keep = sum(releases_per_year) - vintages_to_exclude;

    # Initialise forecasts and outturn
    point_forecasts = Array{Float64}(undef, h, n, releases_to_keep);
    rw_forecasts    = Array{Float64}(undef, h, n, releases_to_keep);
    outturn         = Array{Float64}(undef, h, n, releases_to_keep);

    # Initialise output states
    output_gap       = Array{Float64}(undef, T, releases_to_keep);
    potential_output = Array{Float64}(undef, T, releases_to_keep);
    monthly_gdp      = Array{Float64}(undef, T, releases_to_keep);

    # Initialise remaining states
    BC_clean = Array{Float64}(undef, T, releases_to_keep);
    EP_clean = Array{Float64}(undef, T, releases_to_keep);
    BC       = Array{Float64}(undef, T, releases_to_keep);
    EP       = Array{Float64}(undef, T, releases_to_keep);
    T_INFL   = Array{Float64}(undef, T, releases_to_keep);

    # Loop over the chunks
    for i=1:nyears

        # Print loading status
        @info("Loading chunk $i out of $nyears");

        # Load current chunk
        raw_results = jldopen("$(model_folder)/results/res_chunk$(i).jld");

        # Current indeces
        if i == 1
            start_ind_i = 1;
        else
            start_ind_i = 1 + cumsum(releases_per_year[1:i-1])[end];
        end
        end_ind_i = cumsum(releases_per_year[1:i])[end];

        if start_ind_i > releases_to_keep
            break;
        elseif end_ind_i > releases_to_keep
            end_ind_i = releases_to_keep;
        end

        # i-th ranges
        ith_vintage_id_range = start_ind_i:end_ind_i;
        ith_input_range = 1:end_ind_i-start_ind_i+1;

        # Store data from current chunk (forecasts)
        point_forecasts[:, :, ith_vintage_id_range] = read(raw_results["point_forecasts"])[:, :, ith_input_range];
        rw_forecasts[:, :, ith_vintage_id_range]    = read(raw_results["rw_forecasts"])[:, :, ith_input_range];
        outturn[:, :, ith_vintage_id_range]         = read(raw_results["outturn"])[:, :, ith_input_range];

        # Store data from current chunk (output states)
        output_gap[:, ith_vintage_id_range]       = dropdims_median(read(raw_results["output_gap"]))[1:T, ith_input_range];
        potential_output[:, ith_vintage_id_range] = dropdims_median(read(raw_results["potential_output"]))[1:T, ith_input_range];
        monthly_gdp[:, ith_vintage_id_range]      = dropdims_median(read(raw_results["potential_output"]) .* (read(raw_results["output_gap"])/100 .+ 1))[1:T, ith_input_range];

        # Store data from current chunk (remaining states)
        BC_clean[:, ith_vintage_id_range] = dropdims_median(read(raw_results["BC_clean"]))[1:T, ith_input_range];
        EP_clean[:, ith_vintage_id_range] = dropdims_median(read(raw_results["EP_clean"]))[1:T, ith_input_range];
        BC[:, ith_vintage_id_range]       = dropdims_median(read(raw_results["BC"]))[1:T, ith_input_range];
        EP[:, ith_vintage_id_range]       = dropdims_median(read(raw_results["EP"]))[1:T, ith_input_range];
        T_INFL[:, ith_vintage_id_range]   = dropdims_median(read(raw_results["T_INFL"]))[1:T, ith_input_range];

        # Optional routine
        if remove_forecast_path

            # Loop over the vintages
            for v=ith_vintage_id_range

                # If not all nans
                if sum(isnan.(output_gap[:, v])) != T

                    # End of the in-sample period
                    vth_nans = findall(isnan.(output_gap[:, v]));
                    if length(vth_nans) > 0
                        end_iis = vth_nans[1]-h-1;
                    else
                        end_iis = T-h;
                    end

                    # Replace forecast period with nans (output states)
                    output_gap[end_iis:end_iis+h, v]       .= NaN;
                    potential_output[end_iis:end_iis+h, v] .= NaN;
                    monthly_gdp[end_iis:end_iis+h, v]      .= NaN;

                    # Replace forecast period with nans (remaining states)
                    BC_clean[end_iis:end_iis+h, v] .= NaN;
                    EP_clean[end_iis:end_iis+h, v] .= NaN;
                    BC[end_iis:end_iis+h, v]       .= NaN;
                    EP[end_iis:end_iis+h, v]       .= NaN;
                    T_INFL[end_iis:end_iis+h, v]   .= NaN;
                end
            end
        end
    end

    # Create the vintages
    data_vintages     = zeros(size(point_forecasts));
    raw_data_vintages = read(chunk0["data_vintages"]);
    last_data_vintage = raw_data_vintages[releases_to_keep];

    # Loop over the vintages
    for v=1:releases_to_keep

        # Loop over the variables
        for i=1:n

            # v-th vintage, i-th variable
            raw_data_vintage_i_v = raw_data_vintages[v][:,data_order[i]];
            start_ind_i = maximum(findall(.~(ismissing.(raw_data_vintage_i_v))));

            # Raw_data
            last_data_vintage_i_v = last_data_vintage[start_ind_i+1:start_ind_i+h, data_order[i]];
            last_data_vintage_i_v[ismissing.(last_data_vintage_i_v)] .= NaN;
            last_data_vintage_i_v = last_data_vintage_i_v |> Array{Float64, 1};

            # Data vintages (this step orders data_vintages implicitely following data_order)
            data_vintages[:, i, v] .= last_data_vintage_i_v;
        end
    end

    date = sort(unique(read(chunk0["df_vintages"])[!, :date]))[1:T-h];

    @info("Loading complete");

    return point_forecasts, rw_forecasts, outturn, data_vintages, date, h, n, chunk0,
           monthly_gdp, output_gap, potential_output, BC_clean, EP_clean, BC, EP, T_INFL;
end

"""
    se_hz(data_vintages, forecasts)

Compute forecast error.
"""
function se_hz(data_vintages, forecasts)

    h, n = size(forecasts[:, :, 1]);

    # Squared error (per release)
    se = (data_vintages .- forecasts).^2

    # Squared error (per horizon)
    se_hz_array = zeros(h, n);
    for hz=1:h
        se_hz_array[hz, :] .= NaN;
        for i=1:n
            data_i_hz = se[hz, i, .~isnan.(se[hz, i, :])];
            if length(data_i_hz) > 0
                se_hz_array[hz,i] = mean(data_i_hz);
            end
        end
    end

    return se_hz_array;
end
