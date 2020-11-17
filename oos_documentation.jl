# ----------------------------------------------------------------------------------------------------------------------
# Manual user input
# ----------------------------------------------------------------------------------------------------------------------

# File path
model_folder = "./models/baseline_oos/";

# true for baseline, false for restricted
is_baseline = true;


# ----------------------------------------------------------------------------------------------------------------------
# Extract output from OOS chunks
# ----------------------------------------------------------------------------------------------------------------------

# Load deps
include("./oos_documentation_deps.jl")

# Load out-of-sample reconstruction output
point_forecasts, rw_forecasts, outturn, data_vintages, date, h, n, chunk0,
       monthly_gdp, output_gap, potential_output, BC_clean, EP_clean, BC, EP, T_INFL = load_oos_recon(16, model_folder, is_baseline);

# Generate vector of unique releases
unique_releases = sort(unique(chunk0["df_vintages"][!, :vintage_id]));

# Compute the SE arrays
TC_SE = se_hz(data_vintages, point_forecasts);
RW_SE = se_hz(data_vintages, rw_forecasts);


# ----------------------------------------------------------------------------------------------------------------------
# MSFE chart and csv output
# ----------------------------------------------------------------------------------------------------------------------

# Initialise chart
titles = chunk0["MNEMONIC"];
splots = Array{Any,1}(undef, length(titles));
gr()

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


# ----------------------------------------------------------------------------------------------------------------------
# Output gap chart and csv output
# ----------------------------------------------------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------------------------------------------------
# Potential output (YoY%) chart and csv output
# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------
# Inflation forecast chart and csv output
# ----------------------------------------------------------------------------------------------------------------------

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
