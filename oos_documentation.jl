# ----------------------------------------------------------------------------------------------------------------------
# Manual user input
# ----------------------------------------------------------------------------------------------------------------------

# Number of (subsequent) oos years to evaluate
nyears = 16;

# File path
model_folder = "./models/baseline_oos";

# true for baseline, false for restricted
is_baseline = true;

# Remove forecast path in the states
remove_forecast_path = true;

# Number of vintages to exclude from the output
vintages_to_exclude = 0;

# File name
output_file_name = "baseline";


# ----------------------------------------------------------------------------------------------------------------------
# Extract output from OOS chunks
# ----------------------------------------------------------------------------------------------------------------------

# Load deps
include("./oos_documentation_deps.jl");

# Load out-of-sample reconstruction output
point_forecasts, rw_forecasts, outturn, data_vintages, date, h, n, chunk0,
       monthly_gdp, output_gap, potential_output, BC_clean, EP_clean, BC, EP, T_INFL = load_oos_recon(nyears, model_folder, is_baseline, remove_forecast_path=remove_forecast_path, vintages_to_exclude=vintages_to_exclude);

# Generate vector of unique releases
unique_releases = sort(unique(read(chunk0["df_vintages"])[!, :vintage_id]))[1:end-vintages_to_exclude];

@info("Creating charts and csv output.");


# ----------------------------------------------------------------------------------------------------------------------
# MSFE csv output
# ----------------------------------------------------------------------------------------------------------------------

# Compute the SE arrays
TC_SE = se_hz(data_vintages, point_forecasts);
RW_SE = se_hz(data_vintages, rw_forecasts);

# Save to csv
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_msfe_tc.csv", DataFrame(TC_SE));
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_msfe_rw.csv", DataFrame(RW_SE));


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

    push!(date_ext, Date("01/$(last_month)/$(last_year)", "dd/mm/yyyy"));
end

traces = Array{Any}(undef, size(output_gap, 2));

for i=1:length(traces)
    traces[i] = scatter(x=date_ext, y=output_gap[:,i], showlegend=false);
end

traces = vcat(traces...);

layout = Layout(title="Monthly output gap", titlefont_size=16,
                xaxis=attr(tickfont_size=10, showgrid=true, linecolor="black", mirror=true, nticks=10, tickangle=0, titlefont=attr(size=10)),
                yaxis=attr(zeroline=true, tickfont_size=10, showgrid=true, linecolor="black", nticks=10, mirror=true, range=[-16, 6], titlefont=attr(size=10), title="Percent"));

fig = plot(traces, layout);

# Size
fig.plot.layout["width"]  = 1000;
fig.plot.layout["height"] = 400;

# Margins
fig.plot.layout["margin"][:b]  = 40;
fig.plot.layout["margin"][:t]  = 40;
fig.plot.layout["margin"][:r]  = 40;
fig.plot.layout["margin"][:l]  = 40;

fig.plot.layout["legend"] = attr(y=-0.1, x=0.415, font=attr(size=10), orientation="h");

PlotlyBase.savefig(fig, "$(model_folder)/img/$(output_file_name)_output_gap.pdf", format="pdf");

# Save to csv
df_output_gap_data = [date_ext output_gap];
df_output_gap_names = vcat(:ref_period, Symbol.(unique_releases));
df_output_gap = DataFrame(df_output_gap_data);
rename!(df_output_gap, df_output_gap_names);
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_output_gap.csv", df_output_gap);


# ----------------------------------------------------------------------------------------------------------------------
# Potential output (YoY%) chart and csv output
# ----------------------------------------------------------------------------------------------------------------------

potential_output_yoy = 100*((potential_output[13:end,:] ./ potential_output[1:end-12,:]) .- 1);

traces = Array{Any}(undef, size(potential_output_yoy, 2));

for i=1:length(traces)
    traces[i] = scatter(x=date_ext[13:end], y=potential_output_yoy[:,i], showlegend=false);
end

traces = vcat(traces...);

layout = Layout(title="Monthly potential output (YoY, %)", titlefont_size=16,
                xaxis=attr(tickfont_size=10, showgrid=true, linecolor="black", mirror=true, nticks=10, tickangle=0, titlefont=attr(size=10)),
                yaxis=attr(zeroline=false, tickfont_size=10, showgrid=true, linecolor="black", nticks=10, mirror=true, range=[-2, 6], titlefont=attr(size=10), title="Percent"));

fig = plot(traces, layout);

# Size
fig.plot.layout["width"]  = 1000;
fig.plot.layout["height"] = 400;

# Margins
fig.plot.layout["margin"][:b]  = 40;
fig.plot.layout["margin"][:t]  = 40;
fig.plot.layout["margin"][:r]  = 40;
fig.plot.layout["margin"][:l]  = 40;

fig.plot.layout["legend"] = attr(y=-0.1, x=0.415, font=attr(size=10), orientation="h");

PlotlyBase.savefig(fig, "$(model_folder)/img/$(output_file_name)_potential_output_yoy.pdf", format="pdf");

# Save to csv
df_potential_output_data = [date_ext[13:end] potential_output_yoy];
df_potential_output_names = vcat(:ref_period, Symbol.(unique_releases));
df_potential_output = DataFrame(df_potential_output_data);
rename!(df_potential_output, df_potential_output_names);
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_potential_output_yoy.csv", df_potential_output);


# ----------------------------------------------------------------------------------------------------------------------
# Inflation forecast chart and csv output
# ----------------------------------------------------------------------------------------------------------------------

df_headline_fc = DataFrame([collect(1:36) point_forecasts[:, 6 + ~is_baseline, :]]);
rename!(df_headline_fc, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_forecast_headline.csv", df_headline_fc);

df_headline_outturn = DataFrame([collect(1:36) outturn[:, 6 + ~is_baseline, :]]);
rename!(df_headline_outturn, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_outturn_headline.csv", df_headline_outturn);

for hz=1:h

    # Data (including nans)
    y_hz = df_headline_outturn[hz, 2:end] |> Array;
    fc_hz = df_headline_fc[hz, 2:end] |> Array;

    # Find last observed measurement
    hzth_nans = findall(isnan.(y_hz));
    if length(hzth_nans) > 0
        last_obs_hz  = hzth_nans[1]-1;
    else
        last_obs_hz = size(y_hz, 1);
    end

    # Data (excluding nans)
    y_hz  = y_hz[1:last_obs_hz];
    fc_hz = fc_hz[1:last_obs_hz];

    if hz == 1
        label_name = "$hz month ahead forecast";
    else
        label_name = "$hz months ahead forecast";
    end

    trace1 = scatter(x=unique_releases[1:last_obs_hz], y=fc_hz, line=attr(width=1.4, color=c1), mode="lines", showlegend=true, name=label_name);
    trace2 = scatter(x=unique_releases[1:last_obs_hz], y=y_hz, line=attr(width=1.4, color="black", dash="dot"), mode="lines", showlegend=true, name="Outturn");
        
    layout = Layout(title="Headline inflation (YoY, %)", titlefont_size=12,
                    xaxis=attr(tickfont_size=10, showgrid=true, linecolor="black", mirror=true, nticks=10, tickangle=0, range=[unique_releases[1], unique_releases[last_obs_hz]], titlefont=attr(size=10), title="Releases"),
                    yaxis=attr(zeroline=false, tickfont_size=10, showgrid=true, linecolor="black", nticks=10, mirror=true, titlefont=attr(size=10), title="Percent"));
            
    fig = plot([trace1, trace2], layout);

    # Size
    fig.plot.layout["width"]  = 1000;
    fig.plot.layout["height"] = 400;

    # Margins
    fig.plot.layout["margin"][:b]  = 40;
    fig.plot.layout["margin"][:t]  = 40;
    fig.plot.layout["margin"][:r]  = 40;
    fig.plot.layout["margin"][:l]  = 40;

    fig.plot.layout["legend"] = attr(y=-0.2, x=0.36, font=attr(size=10), orientation="h")

    PlotlyBase.savefig(fig, "$(model_folder)/img/$(output_file_name)_headline_forecast_h$(hz).pdf", format="pdf");
end


# ----------------------------------------------------------------------------------------------------------------------
# Real GDP forecast
# ----------------------------------------------------------------------------------------------------------------------

df_gdp_fc = DataFrame([collect(1:36) point_forecasts[:, 1 + ~is_baseline, :]]);
rename!(df_gdp_fc, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_forecast_gdp.csv", df_gdp_fc);

df_gdp_outturn = DataFrame([collect(1:36) outturn[:, 1 + ~is_baseline, :]]);
rename!(df_gdp_outturn, vcat(:fc_horizon, df_output_gap_names[2:end]));
FileIO.save("$(model_folder)/results_csv/$(output_file_name)_outturn_gdp.csv", df_gdp_outturn);

for hz=3:3:h

    # Data (including nans)
    y_hz = df_gdp_outturn[hz, 2:end] |> Array;
    fc_hz = df_gdp_fc[hz, 2:end] |> Array;

    # Find last observed measurement
    hzth_nans = findall(isnan.(y_hz));
    if length(hzth_nans) > 0
        last_obs_hz  = hzth_nans[1]-1;
    else
        last_obs_hz = size(y_hz, 1);
    end

    # Data (excluding nans)
    y_hz  = y_hz[1:last_obs_hz];
    fc_hz = fc_hz[1:last_obs_hz];

    if hz == 3
        label_name = "$(Int64(hz/3)) quarter ahead forecast";
    else
        label_name = "$(Int64(hz/3)) quarters ahead forecast";
    end

    trace1 = scatter(x=unique_releases[1:last_obs_hz], y=fc_hz, line=attr(width=1.4, color=c1), mode="lines", showlegend=true, name=label_name);
    trace2 = scatter(x=unique_releases[1:last_obs_hz], y=y_hz, line=attr(width=1.4, color="black", dash="dot"), mode="lines", showlegend=true, name="Outturn");
        
    layout = Layout(title="Real GDP", titlefont_size=12,
                    xaxis=attr(tickfont_size=10, showgrid=true, linecolor="black", mirror=true, nticks=10, tickangle=0, range=[unique_releases[1], unique_releases[last_obs_hz]], titlefont=attr(size=10), title="Releases"),
                    yaxis=attr(zeroline=false, tickfont_size=10, showgrid=true, linecolor="black", nticks=10, mirror=true, titlefont=attr(size=10), title="Bil. Chn. 2012\$"));
            
    fig = plot([trace1, trace2], layout);

    # Size
    fig.plot.layout["width"]  = 1000;
    fig.plot.layout["height"] = 400;

    # Margins
    fig.plot.layout["margin"][:b]  = 40;
    fig.plot.layout["margin"][:t]  = 40;
    fig.plot.layout["margin"][:r]  = 40;
    fig.plot.layout["margin"][:l]  = 40;

    fig.plot.layout["legend"] = attr(y=-0.2, x=0.36, font=attr(size=10), orientation="h")

    PlotlyBase.savefig(fig, "$(model_folder)/img/$(output_file_name)_gdp_forecast_h$(Int64(hz/3)).pdf", format="pdf");
end


# ----------------------------------------------------------------------------------------------------------------------
# Save forecast summary to disk
# ----------------------------------------------------------------------------------------------------------------------

@info("Save forecast output summary to disk");
save("$(model_folder)/results/forecast_output_summary.jld", Dict("point_forecasts" => point_forecasts, 
     "rw_forecasts" => rw_forecasts, "outturn" => outturn, "date" => date, "unique_releases" => unique_releases, 
     "h" => h));