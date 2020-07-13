"""
    get_fred_vintages(fred_data_path::String, start_sample::Date, end_sample::Date, oos_start_date::Date)

Build data vintages from Fred data.

The code depends on FredData.jl. Details on the parameters can be found at the link below.
https://research.stlouisfed.org/docs/api/fred/series_observations.html#Description
"""
function get_fred_vintages(fred_data_path::String, start_sample::Date, end_sample::Date, oos_start_date::Date)

    # Initialise FredData instance
    f = Fred(); # FredData must be correctly setup to work

    # Load from Excel data
    excel_input = DataFrame(XLSX.readtable(fred_data_path, "alfred")...);
    tickers = excel_input[!,:fred_tickers] |> Array{String,1};
    MNEMONIC = excel_input[!,:MNEMONIC] |> Array{String,1};
    frequency = excel_input[!,:frequency] |> Array{String,1};
    unit = excel_input[!,:unit] |> Array{String,1};
    aggregation = excel_input[!,:aggregation] |> Array{String,1};

    # Initialise final output
    fred_vintages = DataFrame();

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

    # Return output
    return fred_vintages;
end

"""
    get_local_vintages(local_data_path::String, oos_start_date::Date)

Build data vintages from local data.
"""
function get_local_vintages(local_data_path::String, oos_start_date::Date)

    # Initialise final output
    local_vintages = DataFrame();

    # Load local data
    data, date, nM_local, nQ_local, MNEMONIC_local = read_data(local_data_path);

    # Load local calendar of releases
    raw_calendar = DataFrame(XLSX.readtable(local_data_path, "calendar")...);
    reference    = vcat([Dates.Date(names(raw_calendar)[1+i]) for i=1:size(raw_calendar,2)-1]...);
    releases     = raw_calendar[:, 2:end] |> Array{Union{Missing,Date},2};

    # Match reference with date
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

    # Cut reference and releases
    reference = reference[ind_dates];
    releases = releases[:, ind_dates];

    # Initialise final output
    local_vintages = DataFrame();

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

    # Return output
    return local_vintages;
end

"""
    get_vintages(df_vintages::DataFrame, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{String,1}, h::Int64)

Construct vintages in Array format and transform the data.
"""
function get_vintages(df_vintages::DataFrame, MNEMONIC::Array{String,1}, transf::Array{Int64,1}, transf_annex::Array{String,1}, h::Int64);

    # Re-order the data so that the quarterly series follow the monthly ones
    df_vintages = df_vintages[!, vcat([:vintage_id, :date], Symbol.(MNEMONIC))];
    sort!(df_vintages, [:vintage_id]);

    # Convert to dataflow
    unique_releases = unique(df_vintages[!,:vintage_id]);
    unique_reference = sort(unique(df_vintages[!,:date]));
    data_vintages = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);

    # Loop over the releases
    for i=1:length(unique_releases)

        # Positional index
        ith_position = findall(df_vintages[!,:vintage_id] .== unique_releases[i]);

        # Current vintage (raw version)
        current_vintage_raw = df_vintages[ith_position, 2:end];
        sort!(current_vintage_raw, [:date]);

        # Update vintage
        if i > 1

            # Latest vintage
            latest_vintage = DataFrame(data_vintages[i-1])
            latest_vintage = hcat(DataFrame(:date => unique_reference[1:size(data_vintages[i-1], 1)]), latest_vintage);
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
            sort!(current_vintage, [:date]);

            # Store current_vintage
            data_vintages[i] = [current_vintage[:, 2:end]; missing.*ones(h, length(MNEMONIC))];

        else
            data_vintages[i] = [current_vintage_raw[:, 2:end]; ; missing.*ones(h, length(MNEMONIC))];
        end
    end

    # Return output
    return data_vintages;
end
