# ----------------------------------------------------------------------------------------------------------------------
# Initial settings
# ----------------------------------------------------------------------------------------------------------------------

using JLD, FileIO;
using DataFrames;
using XLSX;
using Plots, Measures;
using Dates, Statistics;


# ----------------------------------------------------------------------------------------------------------------------
# Custom dependencies
# ----------------------------------------------------------------------------------------------------------------------

"""
    dropdims_median(X::Array{Float64,3})

Take median and remove unit dimension.
"""
dropdims_median(X::Array{Float64,3}) = dropdims(median(X, dims=2), dims=2);

"""
    load_data(nyears, model_folder, is_baseline)

Load out-of-sample reconstruction results from the output chunks.
"""
function load_data(nyears, model_folder, is_baseline)

    # Load chunk0
    chunk0 = load("$(model_folder)/results/res_chunk0.jld");

    # Initial settings
    data_order        = chunk0["data_order"];
    releases_per_year = chunk0["releases_per_year"];

    # Dimensions
    h    = size(read(JLD.jldopen("$(model_folder)/results/res_chunk1.jld")["point_forecasts"]))[1];
    T, n = size(chunk0["data"]);

    # Initialise forecasts and outturn
    point_forecasts = Array{Float64}(undef, h, n, sum(releases_per_year));
    rw_forecasts    = Array{Float64}(undef, h, n, sum(releases_per_year));
    outturn         = Array{Float64}(undef, h, n, sum(releases_per_year));

    # Initialise output states
    output_gap       = Array{Float64}(undef, T, sum(releases_per_year));
    potential_output = Array{Float64}(undef, T, sum(releases_per_year));
    monthly_gdp      = Array{Float64}(undef, T, sum(releases_per_year));

    # Initialise remaining states
    BC_clean = Array{Float64}(undef, T, sum(releases_per_year));
    EP_clean = Array{Float64}(undef, T, sum(releases_per_year));
    BC       = Array{Float64}(undef, T, sum(releases_per_year));
    EP       = Array{Float64}(undef, T, sum(releases_per_year));
    T_INFL   = Array{Float64}(undef, T, sum(releases_per_year));

    # Loop over the chunks
    for i=1:nyears

        # Load current chunk
        raw_results = load("$(model_folder)/results/res_chunk$(i).jld");

        # Current indeces
        if i == 1
            start_ind_i = 1;
        else
            start_ind_i = 1 + cumsum(releases_per_year[1:i-1])[end];
        end
        end_ind_i = cumsum(releases_per_year[1:i])[end];

        # Store data from current chunk (forecasts)
        point_forecasts[:,:,start_ind_i:end_ind_i] = raw_results["point_forecasts"];
        rw_forecasts[:,:,start_ind_i:end_ind_i]    = raw_results["rw_forecasts"];
        outturn[:,:,start_ind_i:end_ind_i]         = raw_results["outturn"];

        # Store data from current chunk (output states)
        output_gap[:,start_ind_i:end_ind_i]       = dropdims_median(raw_results["output_gap"]);
        potential_output[:,start_ind_i:end_ind_i] = dropdims_median(raw_results["potential_output"]);
        monthly_gdp[:,start_ind_i:end_ind_i]      = dropdims_median(raw_results["output_gap"] .+ raw_results["potential_output"]);

        # Store data from current chunk (remaining states)
        BC_clean[:,start_ind_i:end_ind_i] = dropdims_median(raw_results["BC_clean"]);
        EP_clean[:,start_ind_i:end_ind_i] = dropdims_median(raw_results["EP_clean"]);
        BC[:,start_ind_i:end_ind_i]       = dropdims_median(raw_results["BC"]);
        EP[:,start_ind_i:end_ind_i]       = dropdims_median(raw_results["EP"]);
        T_INFL[:,start_ind_i:end_ind_i]   = dropdims_median(raw_results["T_INFL"]);
    end

    # Create the vintages
    data_vintages     = zeros(size(point_forecasts));
    raw_data_vintages = chunk0["data_vintages"];
    last_data_vintage = raw_data_vintages[end];

    # Loop over the vintages
    for v=1:size(data_vintages,3)

        # Loop over the variables
        for i=1:n

            # v-th vintage, i-th variable
            raw_data_vintage_i_v = raw_data_vintages[v][:,data_order[i]]
            start_ind_i = maximum(findall(.~(ismissing.(raw_data_vintage_i_v))));

            # Raw_data
            last_data_vintage_i_v = last_data_vintage[start_ind_i+1:start_ind_i+h, data_order[i]];
            last_data_vintage_i_v[ismissing.(last_data_vintage_i_v)] .= NaN;
            last_data_vintage_i_v = last_data_vintage_i_v |> Array{Float64, 1};

            # Data vintages (this step orders data_vintages implicitely following data_order)
            data_vintages[:,i,v] .= last_data_vintage_i_v;
        end
    end

    date = sort(unique(chunk0["df_vintages"][!, :date]));

    return point_forecasts, rw_forecasts, outturn, data_vintages, date, h, n, chunk0,
           monthly_gdp, output_gap, BC_clean, EP_clean, BC, EP, T_INFL;
end

"""
    se_hz(data_vintages, forecasts)

Compute forecast error.
"""
function se_hz(data_vintages, forecasts)

    h, n = size(forecasts[:,:,1]);

    # Squared error (per release)
    se = (data_vintages .- forecasts).^2

    # Squared error (per horizon)
    se_hz_array = zeros(h, n);
    for hz=1:h
        se_hz_array[hz,:] .= NaN;
        for i=1:n
            data_i_hz = se[hz,i,.~isnan.(se[hz,i,:])];
            if length(data_i_hz) > 0
                se_hz_array[hz,i] = mean(data_i_hz);
            end
        end
    end

    return se_hz_array;
end
