using JLD, FileIO;
using DataFrames;
using XLSX;
using Plots, Measures;
using Dates, Statistics;

function load_data(nyears)

    # Load Results
    chunk0 = load("./results/res_chunk0.jld");
    data_order = chunk0["data_order"];

    point_forecasts=Array{};
    outturn=Array{};
    rw_forecasts=Array{};
    monthly_gdp=Array{};
    output_gap=Array{};
    BC_clean=Array{};
    EP_clean=Array{};
    BC=Array{};
    EP=Array{};
    T_INFL=Array{};

    for i=1:nyears
        raw_results = load("./results/res_chunk$(i).jld");

        if i == 1
            point_forecasts=raw_results["point_forecasts"];
            outturn=raw_results["outturn"];
            rw_forecasts=raw_results["rw_forecasts"];

            monthly_gdp = dropdims(median(raw_results["monthly_gdp"], dims=2), dims=2);
            output_gap  = dropdims(median(raw_results["output_gap"], dims=2), dims=2);
            BC_clean    = dropdims(median(raw_results["BC_clean"], dims=2), dims=2);
            EP_clean    = dropdims(median(raw_results["EP_clean"], dims=2), dims=2);
            BC          = dropdims(median(raw_results["BC"], dims=2), dims=2);
            EP          = dropdims(median(raw_results["EP"], dims=2), dims=2);
            T_INFL      = dropdims(median(raw_results["T_INFL"], dims=2), dims=2);

        else
            point_forecasts=cat(dims=3,point_forecasts,raw_results["point_forecasts"]);
            outturn=cat(dims=3,outturn,raw_results["outturn"]);
            rw_forecasts=cat(dims=3,rw_forecasts,raw_results["rw_forecasts"]);

            monthly_gdp = cat(dims=2, monthly_gdp, dropdims(median(raw_results["monthly_gdp"], dims=2), dims=2));
            output_gap  = cat(dims=2, output_gap, dropdims(median(raw_results["output_gap"], dims=2), dims=2));
            BC_clean    = cat(dims=2, BC_clean, dropdims(median(raw_results["BC_clean"], dims=2), dims=2));
            EP_clean    = cat(dims=2, EP_clean, dropdims(median(raw_results["EP_clean"], dims=2), dims=2));
            BC          = cat(dims=2, BC, dropdims(median(raw_results["BC"], dims=2), dims=2));
            EP          = cat(dims=2, EP, dropdims(median(raw_results["EP"], dims=2), dims=2));
            T_INFL      = cat(dims=2, T_INFL, dropdims(median(raw_results["T_INFL"], dims=2), dims=2));
        end
    end

    # create the vintages
    data_vintages = zeros(size(point_forecasts));
    raw_data_vintages = chunk0["data_vintages"];
    last_data_vintage = raw_data_vintages[end];

    h, n = size(point_forecasts[:,:,1]);

    # Loop over the  vintages
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

# Run documentation
point_forecasts, rw_forecasts, outturn, data_vintages, date, h, n, chunk0,
       monthly_gdp, output_gap, BC_clean, EP_clean, BC, EP, T_INFL = load_data(16);

# Generate unique releases from DataFrame
unique_releases = sort(unique(chunk0["df_vintages"][!, :vintage_id]));

# Initialise SE arrays
TC_SE = se_hz(data_vintages, point_forecasts);
RW_SE = se_hz(data_vintages, rw_forecasts);

# Create charts
titles = chunk0["MNEMONIC"];
splots = Array{Any,1}(undef, length(titles));
gr()

"""
MSFE
"""

for i=1:length(titles)

    ind_not_nan = .~isnan.(TC_SE[:,i]);
    if sum(ind_not_nan) != length(ind_not_nan)
        xaxis_name = "Forecast horizon (quarters)";
        xaxis_label = 3:3:12;
    else
        xaxis_name = "Forecast horizon (months)";
        xaxis_label = 3:3:36;
    end

    if i != length(titles)
        splots[i] = plot(TC_SE[ind_not_nan,i], title=titles[i], label="TC model", legend=false, color="red", framestyle=:box, titlefont=font(10), xaxis=(xaxis_name, xaxis_label, font(8)), yaxis=("MSE", font(8)));
        plot!(RW_SE[ind_not_nan,i], label="RW benchmark", legend=false, color="green");
    else
        splots[i] = plot(TC_SE[ind_not_nan,i], title=titles[i], label="TC model", legend=:bottomright, color="red", framestyle=:box, titlefont=font(10), xaxis=(xaxis_name, xaxis_label, font(8)), yaxis=("MSE", font(8)));
        plot!(RW_SE[ind_not_nan,i], label="RW benchmark", legend=:bottomright, color="green");
    end
end

p1 = plot(splots..., layout=(4,2), size=(1200,1000));
Plots.savefig(p1, "./img/rmsfe_plot.pdf");

# Save to csv
FileIO.save("./results_csv/msfe_tc.csv", DataFrame(TC_SE));
FileIO.save("./results_csv/msfe_rw.csv", DataFrame(RW_SE));


"""
Output gap
"""

# Extend date vector
date_ext = copy(date) |> Array{Date,1};

for hz=1:size(output_gap,1)-length(date)
    last_month = month(date_ext[end]);
    last_year = year(date_ext[end]);

    if last_month == 12
        last_month = 1;
        last_year += 1;
    else
        last_month += 1;
    end

    push!(date_ext, Date("01/$(last_month)/$(last_year)", "dd/mm/yyyy"))
end

p2 = plot(date_ext, output_gap, title="Output gap", legend=false, framestyle=:box, titlefont=font(10), xaxis=(font(8)), yaxis=("Percent", font(8)), size=(600,250));
Plots.savefig(p2, "./img/output_gap.pdf");

# Save to csv
df_output_gap_data = [date_ext output_gap];
df_output_gap_names = vcat(:ref_period, Symbol.(unique_releases));
df_output_gap = DataFrame(df_output_gap_data);
rename!(df_output_gap, df_output_gap_names);
FileIO.save("./results_csv/output_gap.csv", df_output_gap);


"""
Inflation forecast
"""

df_headline_fc = DataFrame([collect(1:36) point_forecasts[:,6,:]]);
rename!(df_headline_fc, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("./results_csv/forecast_headline.csv", df_headline_fc);

df_headline_outturn = DataFrame([collect(1:36) outturn[:,6,:]]);
rename!(df_headline_outturn, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("./results_csv/outturn_headline.csv", df_headline_outturn);

for hz=1:h
    y_hz = df_headline_outturn[hz,2:end] |> Array;
    fc_hz = df_headline_fc[hz,2:end] |> Array;
    p_hz = plot(unique_releases, y_hz, color=:black, label="Headline inflation", framestyle=:box, titlefont=font(10), xaxis=(font(8)), yaxis=("Percent", font(8)), size=(600,250));
    plot!(unique_releases, fc_hz, color=:blue, label="$(hz)-step ahead forecast");
    Plots.savefig(p_hz, "./img/headline_forecast_h$(hz).pdf");
end
