"""
    get_fred_vintages(fred_data_path::String, start_sample::Date, end_sample::Date, oos_start_date::Date)

Build data vintages from Fred data.

The code depends on FredData.jl. Details on the parameters can be found at the link below.
https://research.stlouisfed.org/docs/api/fred/series_observations.html#Description
"""
function get_fred_vintages(fred_data_path::String, start_sample::Date, end_sample::Date, oos_start_date::Date)

    # Initialise final output
    fred_vintages = DataFrame();

    # Initialise FredData instance
    f = Fred(); # FredData must be correctly setup to work

    try

        # Load from Excel data
        excel_input = DataFrame(XLSX.readtable(fred_data_path, "alfred")...);
        tickers = excel_input[!,:fred_tickers] |> Array{String,1};
        MNEMONIC = excel_input[!,:MNEMONIC] |> Array{String,1};
        frequency = excel_input[!,:frequency] |> Array{String,1};
        unit = excel_input[!,:unit] |> Array{String,1};
        aggregation = excel_input[!,:aggregation] |> Array{String,1};

        # Loop over tickers
        for i=1:length(tickers)

            # Slow down the requests - otherwise, FredData might crash
            if i > 1
                sleep(rand(1:0.5:5));
            end

            # Get data from Fred
            ith_data = get_data(f, tickers[i];
                                realtime_start = Dates.format(oos_start_date, "yyyy-mm-dd"),
                                observation_start = Dates.format(start_sample, "yyyy-mm-dd"),
                                observation_end = Dates.format(end_sample, "yyyy-mm-dd"),
                                frequency = frequency[i],
                                units = unit[i]).data[:,[:realtime_start, :date, :value]];

            # Convert NaN to missing
            ith_data[!,:value] = convert(Array{Union{Missing, Float64}}, ith_data[!,:value]);
            ith_data[isnan.(ith_data[!,:value]),:value] .= missing;

            # Rename value to tickers[i]
            rename!(ith_data, (:value => Symbol(MNEMONIC[i])));
            rename!(ith_data, (:realtime_start => :vintage_id));

            #=
            Quarterly data is not aligned correctly at source (e.g., Q1 has Jan as reference period).
            This block of code fixes it.
            =#
            if frequency[i] == "q"
                ith_data[!,:date] = Dates.firstdayofmonth.(ith_data[!,:date]);
                ith_data[!,:date] .+= Month(2);
            end

            # Set reference periods to eom format
            ith_data[!,:date] = Dates.lastdayofmonth.(ith_data[!,:date]);

            # Store new data
            if i == 1
                fred_vintages = copy(ith_data);
            else
                fred_vintages = outerjoin(fred_vintages, ith_data, on=[:date, :vintage_id]);
            end
        end

    catch
        @warn("get_fred_vintages crashed. \\ Is $fred_data_path empty? If not, double check the input file.\\ df_fred_vintages is defined as an empty DataFrame");
    end

    # Chronological order
    sort!(fred_vintages, :vintage_id);

    # Return output
    return fred_vintages;
end

"""
    get_local_vintages(local_data_path::String, end_sample::Date, oos_start_date::Date)

Build data vintages from local data.
"""
function get_local_vintages(local_data_path::String, end_sample::Date, oos_start_date::Date)

    # Initialise final output
    local_vintages = DataFrame();

    try

        # Load local data
        data, date, nM_local, nQ_local, MNEMONIC_local = read_data(local_data_path);

        # Load local calendar of releases
        raw_calendar = DataFrame(XLSX.readtable(local_data_path, "calendar")...);
        reference    = vcat([Dates.Date(names(raw_calendar)[1+i]) for i=1:size(raw_calendar,2)-1]...);
        releases     = raw_calendar[:, 2:end] |> Array{Union{Missing,Date},2};

        # Resize calendar wrt end_sample
        reference = reference[reference .<= end_sample];
        releases = releases[:, reference .<= end_sample];

        # This block removes the calendar entries that do not match any date in `date`
        ind_dates = [];
        ind_dataset = [];
        for i=1:length(reference)
            ith_position = findall(date .== reference[i]);
            if length(ith_position) > 0
                push!(ind_dates, i);
                push!(ind_dataset, ith_position[1]);
            end
        end

        min_ind_dates = minimum(ind_dates);
        min_ind_dataset = minimum(ind_dataset);

        reference = reference[ind_dates];
        releases = releases[:, ind_dates];

        # Loop over local series
        for i=1:nM_local+nQ_local

            # Data excluding initial datapoints
            ith_values = Array{Union{Missing, Float64}}(data[ind_dataset,i]);
            ith_data = DataFrame(:vintage_id => releases[i,:], :date => reference, Symbol(MNEMONIC_local[i]) => ith_values);

            # Initial datapoints (pre-calendar)
            if min_ind_dataset > 1
                ith_values_init = Array{Union{Missing, Float64}}(data[1:min_ind_dataset-1,i]);
                ith_releases_init = [minimum(skipmissing(ith_data[!,:vintage_id])) for j=1:min_ind_dataset-1];
                ith_data_init =  DataFrame(:vintage_id => ith_releases_init, :date => date[1:min_ind_dataset-1], Symbol(MNEMONIC_local[i]) => ith_values_init);
                ith_data = vcat(ith_data_init, ith_data);
            end

            # Skip releases of missings (i.e., fake releases for quarterly data)
            ith_data = ith_data[findall(.~ismissing.(ith_data[!,Symbol(MNEMONIC_local[i])])), :];

            # Correct vintage_id
            ith_data[ith_data[!, :vintage_id] .< oos_start_date, :vintage_id] .= oos_start_date;

            # Set reference periods to eom format
            ith_data[!,:date] = Dates.lastdayofmonth.(ith_data[!,:date]);

            # Store new data
            if i == 1
                local_vintages = copy(ith_data);
            else
                local_vintages = outerjoin(local_vintages, ith_data, on=[:date, :vintage_id]);
            end
        end

    catch
        @warn("get_local_vintages crashed. \\ Is $local_data_path empty? If not, double check the input file.\\ df_local_vintages is defined as an empty DataFrame");
    end

    # Chronological order
    sort!(local_vintages, :vintage_id);

    # Return output
    return local_vintages;
end

"""
    transform_vintage(vintage::DataFrame, start_sample::Date, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{Any,1}, nM::Int64)

Transform individual data vintage.
"""
function transform_vintage(vintage::DataFrame, start_sample::Date, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{Any,1}, nM::Int64)

    # Loop over the columns to transform
    for i=findall(transf .!= 0)

        ith_symbol = Symbol(MNEMONIC[i]);

        # YoY % transformation
        if transf[i] == 1
            vintage[13:end, ith_symbol] .= 100*(vintage[13:end, ith_symbol] ./ vintage[1:end-12, ith_symbol] .- 1);
            vintage[1:12, ith_symbol] .= missing;

        # Rebase and de-annualise
        elseif transf[i] == 2

            # Conversion factor
            ith_conversion = transf_annex[i] / vintage[vintage[!,:date].==lastdayofquarter(start_sample), ith_symbol][1]

            # Rabase and de-annualise
            vintage[:, ith_symbol] .*= ith_conversion/4;

        # Compute SPF expectation in levels, starting from the growth rates (wrt the variable specified in TRANSF_ANNEX)
        elseif transf[i] == 3
            if i <= nM
                error("Transformation not supported for monthly data yet. The transform_vintage function must be updated to enable this feature.")
            else

                # Identify first quarter
                first_not_missing_spf = findall(.~ismissing.(vintage[:, ith_symbol]))[1];

                # Compute growth factor (i.e., 1+growth rate)
                vintage[:, ith_symbol] .+= 1;

                # Compute SPF in levels
                vintage[first_not_missing_spf:end, ith_symbol] .*= vintage[first_not_missing_spf-3:end-3, Symbol(transf_annex[i])];
            end

        # Cycle transformation (2 + compute cycle subtracting trend from variable specified in TRANSF_ANNEX)
        elseif transf[i] == 4
            vintage[:, ith_symbol] .-= vintage[:, Symbol(transf_annex[i])];
            vintage[:, ith_symbol] .*= -1;
        end

        if transf[i] == 3 || transf[i] == 4
            if findall(MNEMONIC .== transf_annex[i])[1] > i && transf[findall(MNEMONIC .== transf_annex[i])[1]] == 2
                error("The $i-th variable must be positioned after $(transf_annex[i]). \\ If $i-th is monthly the transform_vintage function must be updated, or the transformation type changed.")
            end
        end
    end

    # If at least one variable is transformed in YoY% increase the sample of one year skip the first 12 observations
    if sum(transf .== 1) > 0
        vintage = vintage[13:end, :];
    elseif sum(transf .== 1) == 0 && sum(transf .== 3) > 0
        vintage = vintage[4:end, :];
    end

    # Return output
    return vintage;
end

"""
    get_vintages(df_vintages::DataFrame, start_sample::Date, end_sample::Date, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{Any,1}, nM::Int64, h::Int64)

Construct vintages in Array format and transform the data.
"""
function get_vintages(df_vintages::DataFrame, start_sample::Date, end_sample::Date, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{Any,1}, nM::Int64, h::Int64)

    # Re-order the data so that the quarterly series follow the monthly ones
    df_vintages = df_vintages[!, vcat([:vintage_id, :date], Symbol.(MNEMONIC))];
    sort!(df_vintages, :vintage_id);

    # Convert to dataflow
    unique_releases = unique(df_vintages[!,:vintage_id]);
    unique_reference = sort(unique(df_vintages[!,:date]));
    data_vintages = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);
    data_vintages_untransformed = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);

    # Date sample range
    date_sample_range = Dates.lastdayofmonth.(collect(start_sample:Dates.Month(1):end_sample));

    # Loop over the releases
    for i=1:length(unique_releases)
        if i == 1 || mod(i, 100) == 0
            println("Build vintage $i (out of $(length(unique_releases)))")
        end

        # Positional index
        ith_position = findall(df_vintages[!,:vintage_id] .== unique_releases[i]);

        # Current vintage (raw version)
        current_vintage_raw = df_vintages[ith_position, 2:end];
        sort!(current_vintage_raw, :date);

        # Update vintage
        if i == 1
            current_vintage = copy(current_vintage_raw);
        else
            # Latest vintage
            latest_vintage = DataFrame(data_vintages_untransformed[i-1]);
            latest_vintage_obs = size(latest_vintage, 1);
            latest_vintage = hcat(DataFrame(:date => unique_reference[1:latest_vintage_obs]), latest_vintage);
            rename!(latest_vintage, vcat(:date, Symbol.(MNEMONIC)));

            # Potential revisions and new releases for already observed reference months
            revisions = semijoin(current_vintage_raw, latest_vintage, on=:date);

            # `revisions` correspondence in `latest_vintage`
            revisions_correspondence = semijoin(latest_vintage, current_vintage_raw, on=:date);

            # Merge the latest two datasets
            for j=2:size(revisions, 2)
                jth_ind_missings = ismissing.(revisions[:,j]);
                revisions[jth_ind_missings, j] .= revisions_correspondence[jth_ind_missings, j];
            end

            # Non-revised data
            non_revised = antijoin(latest_vintage, current_vintage_raw, on=:date);

            # New reference months
            new_entries = antijoin(current_vintage_raw, latest_vintage, on=:date);

            # Define current_vintage
            current_vintage = vcat(non_revised, revisions, new_entries);
            sort!(current_vintage, :date);
        end

        # Store untransformed vintage
        data_vintages_untransformed[i] = current_vintage[:, 2:end];

        # Add missing rows (if necessary)
        ith_date_sample_range = date_sample_range[date_sample_range .<= maximum(current_vintage[!,:date])];
        current_vintage = outerjoin(DataFrame(:date=>ith_date_sample_range), current_vintage, on=:date);
        sort!(current_vintage, :date);

        # Transform vintage
        current_vintage = transform_vintage(current_vintage, start_sample, MNEMONIC, transf, transf_annex, nM);

        # Store current vintage
        data_vintages[i] = current_vintage[:, 2:end];

        # Add missings to account for the forecast horizon
        data_vintages[i] = vcat(data_vintages[i], missing.*ones(h, length(MNEMONIC)));
    end

    # Build unique_years and data_vintages_year
    vintage_years = [Year(unique_releases[i]).value for i=1:size(unique_releases,1)];
    unique_years = unique(vintage_years);
    releases_per_year = vcat([sum(vintage_years .== unique_years[i]) for i=1:length(unique_years)]...);

    data_vintages_year = Array{Array{Array{Union{Missing,Float64},2}, 1}}(undef, length(unique_years));
    for i=1:length(unique_years)
        data_vintages_year[i] = Array{Array{Union{Missing,Float64},2},1}(undef, releases_per_year[i]);
    end

    ind_year = findall([1; diff(vintage_years)] .== 1);
    ind_year = [ind_year [ind_year[2:end] .- 1; length(vintage_years)]];

    for i=1:size(ind_year, 1)
        for j=ind_year[i,1]:ind_year[i,2]
            data_vintages_year[i][j-ind_year[i,1]+1] = data_vintages[j];
        end
    end

    # Return output
    return data_vintages, data_vintages_year, unique_years, releases_per_year;
end
