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
        rename!(ith_data, (:value => Symbol(tickers[i])));
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

# Based on the function get_dataflow(), gets Dataframe and outputs Arrays of vintages
"""
    get_dataflow_df(data::DataFrame, Transform::Array{String,1}, h::Int64, oos_start_date::DateTime, OptionPrintFinal = "Yes")
        This Function convert a DataFrame into Arrays of Arrays
        Inputs:
            data: DataFrame with data. First Column must be vintages dates. Second Column Observation Dates. Then Vars
            Transform: same size of number of vars. None no transformations, PC1 YoY percentage change.
            h: forecast horizon
            oos_start_date: when the out of sample starts
            OptionPrintFinal: if Yes, then last vintage will be printed in an csv file ("Final_Dataset.csv")
        Outputs:
            data_vintages: Array of vintages stored in Array
            data_vintages_year: Array (years) of Array (vintages per year) of vintages stored in Array
            date_vintages: Array of Obs dates -- linked to data_vintages
            date_vintages_year: Array of Array of Obs dates -- linked to data_vintages_year
            unique_years: Years of the realeses
            releases_per_year: number of releases per year
            calendar_vintages: Array of vintages (calendars) stored in Array
"""
function get_dataflow_df(data::DataFrame, Transform::Array{String,1}, h::Int64, oos_start_date::DateTime, OptionPrintFinal = "Yes")
    println("")
    printstyled("Note: First Column must be vintages dates. Second Column Observation Dates. Then Vars.", color = :yellow)
    println("")
    # Size
    n = size(data,2);
    # Starting of the sample (when growth rates are computed, the first year will be dropped)
    Ind_start = 1;
    # Set all releases before the out of sample exercise to the date equal or greater than the start of the out of sample
    data[data[:,1].<oos_start_date,1] = minimum(data[data[:,1].>=oos_start_date,1]);
    # Data
    raw_calendar = data[:,1];
    releases           = data[:,1];
    Obs_date           = data[:,2];
    # Unique releases and years
    unique_releases = sort(unique(skipmissing(releases[:])));
    years           = vcat([Dates.year(unique_releases[i]) for i=1:length(unique_releases)]...);
    unique_years    = sort(unique(years));
    # Number of releases per year
    releases_per_year  = vcat([sum(years .== unique_years[i]) for i=1:length(unique_years)]...);
    # Initialize variables
    data_vintages      = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);
    data_vintages_year = Array{
        Array{
            Array{Union{Missing,Float64},2}
            ,1}}(undef,length(unique_years));
    for i=1:length(unique_years)
        data_vintages_year[i] = Array{Array{Union{Missing,Float64},2},1}(undef,releases_per_year[i]);
    end
    date_vintages      = Array{Array{Any,1}}(undef,length(unique_releases),1);
    date_vintages_year = Array{
        Array{
            Array{Any,1}
            ,1}}(undef,length(unique_years));
    for i=1:length(unique_years)
        date_vintages_year[i] = Array{Array{Any,1},1}(undef,releases_per_year[i]);
    end
    calendar_vintages = Array{Array{Any,2}}(undef,length(unique_releases));
    # Vintage ID
    v = 1;
    for year=1:length(unique_years)
        for v_year=1:releases_per_year[year]
                # Find the release in "releases"
                release = unique_releases[v];
                # Select already released
                data_v_df = data[releases.<=release,:];
                # Final Data Frame
                final_data_v_df = data_v_df[:,2:2];
                final_data_v_df = unique(final_data_v_df,1);
                max_date = maximum(final_data_v_df[:date]);
                min_date = minimum(final_data_v_df[:date]);
                N_periods = (Dates.year(max_date) -  Dates.year(min_date)+1)*12;
                # Generate the dates
                Alldates_now = Dates.lastdayofmonth.([Date(Dates.year(min_date) + floor((i-1)/12), i - 12*floor((i-1)/12), 1) for i = 1:N_periods]);
                # Drop not needed
                Alldates_now = Alldates_now[Alldates_now.>=min_date];
                Alldates_now = Alldates_now[Alldates_now.<=max_date];
                # Update
                final_data_v_df = DataFrame(date=Alldates_now);
                final_data_v_df = sort(final_data_v_df);
                final_forcalendar_v = unique(final_data_v_df,1);
                final_data_v_df = unique(final_data_v_df,1);
                # Check whether >1 releases are present for the same observation date and if so take the latest release
                for j = 3:size(data_v_df,2)
                    thisdataframe = data_v_df[:,[1,2,j]];
                    # drop missings
                    thisdataframe = thisdataframe[ismissing.(thisdataframe[:,3]).==false,:];
                    # Order them by publication dates: first will be most updated
                    thisdataframe = sort(thisdataframe, 1, rev=true);
                    # drop duplicates except first one, being sorted on publication date!
                    # Could have used unique(.), but not sure the algo will stay the same
                    # in next versions
                    #thisdataframe = unique(thisdataframe,2)
                    alldates = unique(thisdataframe[:,2])
                    Index_Ineed = (zeros(size(alldates))) |> Array{Int64,1}
                    for ij = 1:size(alldates,1)
                        # being sorted in descending order (release date) the minimum index is what I am interested in
                        Index_Ineed[ij] = minimum(findall(thisdataframe[:,2] .== alldates[ij]))
                    end
                    # Keep only last updates
                    thisdataframe = thisdataframe[Index_Ineed,:]
                    # Calendar
                    forcalendar = thisdataframe[:,[2,1]];
                    # Drop Vintage date
                    thisdataframe = thisdataframe[:,[2,3]]
                    # Finally Join it
                    final_data_v_df = join(final_data_v_df, thisdataframe, on = :date, kind = :outer)
                    final_data_v_df = sort(final_data_v_df,1)
                    # Add a number to the name
                    if j>3
                        thisnamesnow = names(forcalendar)
                        thisnamesnow[2] = Symbol(String(thisnamesnow[2]) * "_" * string(j-3))
                        names!(forcalendar, thisnamesnow)
                    end
                    final_forcalendar_v = join(final_forcalendar_v, forcalendar, on = :date, kind = :outer)
                    final_forcalendar_v = sort(final_forcalendar_v,1)
                end
                #final_data_v_array[:,1] = DateTime.date2julian(final_data_v_array[:,1])
                for j = 1:size(Transform,1)
                    #println("transf on var: ", j)
                    IndexThisVar = j + 1;
                    ThisSeries = final_data_v_df[:,IndexThisVar];
                    if Transform[j] == "delta"
                        ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                        ThisSeriesnew[:] .= missing;
                        ThisSeriesnew[2:end] = ThisSeries[2:end] - ThisSeries[1:end-1]
                        #Ind_start = 13;
                        #IndexThisVar = j + 1;
                        #ThisSeries = final_data_v_df[:,IndexThisVar];
                        #ThisSeriesnew = ThisSeries[:];
                        #ThisSeriesnew[:] .= missing;
                        #if size(ThisSeries,1)>12
                        #    # Check the difference is indeed of 12 months
                        #    a = unique(Dates.year.(final_data_v_df[13:end,1]) - Dates.year.(final_data_v_df[1:end-12,1]))
                        #    @assert (size(a,1) == 1) & (a[1] == 1)
                        #    # If it blocks here it is becasue Dates are not alligned, so cannot compute PC1
                        #    ThisSeriesnew[13:end] = (ThisSeries[13:end]./ThisSeries[1:end-12] .- 1)*100;
                        #end
                    elseif Transform[j] == "2delta"
                        for j_transf_run = 1:2
                            ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                            ThisSeriesnew[:] .= missing;
                            ThisSeriesnew[2:end] = ThisSeries[2:end] - ThisSeries[1:end-1]
                            ThisSeries = deepcopy(ThisSeriesnew)
                        end
                    elseif Transform[j] == "log"
                        ThisSeriesnew = log.(deepcopy(ThisSeries))
                    elseif Transform[j] == "delta log"
                        ThisSeries = log.(ThisSeries)
                        ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                        ThisSeriesnew[:] .= missing;
                        ThisSeriesnew[2:end] = ThisSeries[2:end] - ThisSeries[1:end-1]
                    elseif Transform[j] == "2delta log"
                        ThisSeries = log.(ThisSeries)
                        for j_transf_run = 1:2
                            ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                            ThisSeriesnew[:] .= missing;
                            ThisSeriesnew[2:end] = ThisSeries[2:end] - ThisSeries[1:end-1]
                            ThisSeries = deepcopy(ThisSeriesnew)
                        end
                    elseif Transform[j] == "delta log quart"
                        ThisSeries = log.(ThisSeries)
                        ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                        ThisSeriesnew[:] .= missing;
                        ThisSeriesnew[4:end] = ThisSeries[4:end] - ThisSeries[1:end-3]
                    elseif Transform[j] == "None"
                        ThisSeriesnew = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                    elseif Transform[j] == "PC1"
                        ThisSeriesnew = ThisSeries[:];
                        ThisSeriesnew[:] .= missing;
                        #global final_data_v_df_f = final_data_v_df
                        #global ThisSeries_f = ThisSeries
                        if size(ThisSeries,1)>12
                            # Check the difference is indeed of 12 months
                            a = unique(Dates.year.(final_data_v_df[13:end,1]) - Dates.year.(final_data_v_df[1:end-12,1]))
                            @assert (size(a,1) == 1) & (a[1] == 1)
                            # If it blocks here it is becasue Dates are not alligned, so cannot compute PC1
                            ThisSeriesnew[13:end] = (ThisSeries[13:end]./ThisSeries[1:end-12] .- 1)*100;
                        end
                    elseif Transform[j] == "Index quart"
                        g_returns = deepcopy(ThisSeries) |> Array{Union{Missing, Float64},1};
                        missing_structure_first = ismissing.(ThisSeries);
                        g_returns[:] .= 1;
                        g_returns[4:end] = ThisSeries[4:end]./ThisSeries[1:end-3];
                        missing_structure = ismissing.(g_returns);
                        g_returns[missing_structure] .= 1;
                        ThisSeriesnew = 100*cumprod(g_returns) |> Array{Union{Missing, Float64},1};
                        ThisSeriesnew[missing_structure_first] .= missing;
                    end
                    final_data_v_df[:,IndexThisVar] = ThisSeriesnew;
                end
                # print final
                if (OptionPrintFinal == "Yes") & (v == size(unique_releases,1))
                    CSV.write("Final_Dataset.csv", final_data_v_df);
                    printstyled("Last vintage data printed in the csv file: Final_Dataset.csv", color=:green)
                end
                final_data_v_array = convert(Array{Any,2},final_data_v_df[Ind_start:end,:]);
                final_forcalendar_v = convert(Array{Union{Missing,DateTime},2}, final_forcalendar_v[Ind_start:end,:]);
                # Add missings to forecasting horizon
                len_v = size(final_data_v_array[:,2:end],1);
                n_vars = size(final_data_v_array[:,2:end],2);
                data_v = missing.*ones(len_v+h, n_vars) |> Array{Union{Missing,Float64},2};
                data_v[1:len_v,:] = final_data_v_array[:,2:end];
                data_vintages[v]                 = data_v;
                date_vintages[v]                 = convert(Array{Union{Missing,DateTime},1},final_data_v_array[:,1]);
                data_vintages_year[year][v_year] = data_v;
                date_vintages_year[year][v_year] = convert(Array{Union{Missing,DateTime},1},final_data_v_array[:,1]);
                calendar_vintages[v]             = final_forcalendar_v;
                # Update vintage ID
                v += 1;
        end
    end
    return data_vintages, data_vintages_year, date_vintages, date_vintages_year, unique_years, releases_per_year, calendar_vintages;
end


# Functions to generate Arrays for target variables
"""
 THIS FUNCTION NEED TO BE CHECKED
    get_dataflow_target_rev(data::DataFrame, Transform::Array{String,1}, h_years::Int64, oos_start_date::DateTime, Final_Size::Int64, OptionPrintFinal = "Yes")
        This Function convert a DataFrame into Arrays of Arrays!
        It differs from get_dataflow_df, as it includes in each vintage also future values for a prediction horizon defined in "h".
        When multiple releases are present for the same observation, the first one that is higher or equal to the current release date is taken.
        Inputs:
            data: DataFrame with data. First Column must be vintages dates. Second Column Observation Dates. Then Vars
            Transform: same size of number of vars. None no transformations, PC1 YoY percentage change.
            h_years: forecast horizon in years
            oos_start_date: when the out of sample starts
            Final_Size: the size of the latest dataset (the one with most updated vintages)
            OptionPrintFinal: if Yes, then last vintage will be printed in an csv file ("Final_Dataset.csv")
        Outputs:
            data_revtr: Array of vintages stored in Array (all same dimension: Final_Size, but filled with missings)
            data_revtr_year: Array (years) of Array (vintages per year) of vintages stored in Array
            date_revtr: Array of Obs dates -- linked to data_revtr (all same dimension: Final_Size, but filled with missings)
            date_revtr_year: Array of Array of Obs dates -- linked to data_revtr_year
            calendar_revtr: Array of vintages (calendars) stored in Array
"""
function get_dataflow_target_rev(data::DataFrame, Transform::Array{String,1}, h_years::Int64, oos_start_date::DateTime, Final_Size::Int64, OptionPrintFinal = "Yes")
    println("")
    printstyled("Note: First Column must be vintages dates. Second Column Observation Dates. Then Vars.", color=:yellow)
    println("")
    # Size
    n = size(data,2);
    # Starting of the sample (when growth rates are computed, the first year will be dropped)
    Ind_start = 1;
    # Set all releases before the out of sample exercise to the date equal or greater than the start of the out of sample
    data[data[:,1].<oos_start_date,1] = minimum(data[data[:,1].>=oos_start_date,1]);
    # Data
    raw_calendar = data[:,1];
    releases           = data[:,1];
    Obs_date           = data[:,2];
    # Unique releases and years
    unique_releases = sort(unique(skipmissing(releases[:])));
    years           = vcat([Dates.year(unique_releases[i]) for i=1:length(unique_releases)]...);
    unique_years    = sort(unique(years));
    # Number of releases per year
    releases_per_year  = vcat([sum(years .== unique_years[i]) for i=1:length(unique_years)]...);
    # Initialize variables
    data_revtr      = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);
    data_revtr_year = Array{
        Array{
            Array{Union{Missing,Float64},2}
            ,1}}(undef,length(unique_years));
    for i=1:length(unique_years)
        data_revtr_year[i] = Array{Array{Union{Missing,Float64},2},1}(undef,releases_per_year[i]);
    end
    date_revtr      = Array{Array{Any,1}}(undef,length(unique_releases),1);
    date_revtr_year = Array{
        Array{
            Array{Any,1}
            ,1}}(undef,length(unique_years));
    for i=1:length(unique_years)
        date_revtr_year[i] = Array{Array{Any,1},1}(undef,releases_per_year[i]);
    end
    calendar_revtr = Array{Array{Any,2}}(undef,length(unique_releases));
    # Vintage ID
    v = 1;
    for year=1:length(unique_years)
        for v_year=1:releases_per_year[year]
                # Find the release in "releases"
                release = unique_releases[v];
                # To the release date, add the forecasting horizon to get the last observation date
                Obs_date_last = Dates.lastdayofmonth(Date(Dates.year(release)+h_years, Dates.month(release), 1));
                # Select observations needed
                data_v_df = data[data[:,2].<Obs_date_last,:];
                # Final Data Frame
                final_data_v_df = data_v_df[:,2:2];
                final_data_v_df = unique(final_data_v_df,1);
                max_date = maximum(final_data_v_df[:date]);
                min_date = minimum(final_data_v_df[:date]);
                N_periods = (Dates.year(max_date) -  Dates.year(min_date)+1)*12;
                # Generate the dates
                Alldates_now = Dates.lastdayofmonth.([Date(Dates.year(min_date) + floor((i-1)/12), i - 12*floor((i-1)/12), 1) for i = 1:N_periods]);
                # Drop not needed
                Alldates_now = Alldates_now[Alldates_now.>=min_date];
                Alldates_now = Alldates_now[Alldates_now.<=max_date];
                # Update
                final_data_v_df = DataFrame(date=Alldates_now);
                final_data_v_df = sort(final_data_v_df);
                final_forcalendar_v = unique(final_data_v_df,1);
                final_data_v_df = unique(final_data_v_df,1);
                # Check whether >1 releases are present for the same observation date and if so take the first release >= relase
                for j = 3:size(data_v_df,2)
                    thisdataframe = data_v_df[:,[1,2,j]];
                    # drop missings
                    thisdataframe = thisdataframe[ismissing.(thisdataframe[:,3]).==false,:];
                    # Order them by publication dates: first will be less updated
                    thisdataframe = sort(thisdataframe, 1, rev=false);
                    # drop duplicates except first one, being sorted on publication date!
                    # Could have used unique(.), but not sure the algo will stay the same
                    # in next versions
                    #thisdataframe = unique(thisdataframe,2)
                    alldates = unique(thisdataframe[:,2])
                    Index_Ineed = (zeros(size(alldates))) |> Array{Int64,1}
                    for ij = 1:size(alldates,1)
                        # if Observation date is less or equal to release date, then take most updated one as long as the release date is less or equal to release
                        if alldates[ij] <= release
                            # Being sorted in ascending order on release date
                            try
                                Index_Ineed[ij] = maximum(findall( (thisdataframe[:,2] .== alldates[ij]) .& ( thisdataframe[:,1] .< Dates.lastdayofmonth(release)) ))
                            catch
                            end
                            # However, if such an observation is not present (release date is higher), then take first release available
                            if Index_Ineed[ij] == 0
                                Index_Ineed[ij] = minimum(findall( (thisdataframe[:,2] .== alldates[ij])))
                            end
                        else
                        # if Observation date is > release date, than take first available such that release date is >= release
                            # Being sorted in ascending order on release date
                            try
                                Index_Ineed[ij] = minimum(findall( (thisdataframe[:,2] .== alldates[ij]) .& (thisdataframe[:,1] .>= release) ))
                            catch
                            end
                            # otherwise take most updated but with release date <= release
                            if Index_Ineed[ij] == 0
                                Index_Ineed[ij] = maximum(findall( (thisdataframe[:,2] .== alldates[ij]) .& (thisdataframe[:,1] .<= release) ))
                            end
                        end
                    end
                    # Keep only updates selected
                    thisdataframe = thisdataframe[Index_Ineed,:]
                    # Calendar
                    forcalendar = thisdataframe[:,[2,1]];
                    # Drop Vintage date
                    thisdataframe = thisdataframe[:,[2,3]]
                    # Add a number to the name
                    if j>3
                        thisnamesnow = names(forcalendar)
                        thisnamesnow[2] = Symbol(String(thisnamesnow[2]) * "_" * string(j-3))
                        names!(forcalendar, thisnamesnow)
                    end
                    # Finally Join it
                    final_data_v_df = join(final_data_v_df, thisdataframe, on = :date, kind = :outer)
                    final_data_v_df = sort(final_data_v_df,1)
                    final_forcalendar_v = join(final_forcalendar_v, forcalendar, on = :date, kind = :outer)
                    final_forcalendar_v = sort(final_forcalendar_v,1)
                end
                #final_data_v_array[:,1] = DateTime.date2julian(final_data_v_array[:,1])
                for j = 1:size(Transform,1)
                    # Transform to yearly growth
                    if Transform[j] == "PC1"
                        IndexThisVar = j + 1;
                        ThisSeries = final_data_v_df[:,IndexThisVar];
                        ThisSeriesnew = ThisSeries[:];
                        ThisSeriesnew[:] .= missing;
                        if size(ThisSeries,1)>12
                            # Check the difference is indeed of 12 months
                            a = unique(Dates.year.(final_data_v_df[13:end,1]) - Dates.year.(final_data_v_df[1:end-12,1]))
                            @assert (size(a,1) == 1) & (a[1] == 1)
                            # If it blocks here it is becasue Dates are not alligned, so cannot compute PC1
                            ThisSeriesnew[13:end] = (ThisSeries[13:end]./ThisSeries[1:end-12] .- 1)*100;
                        end
                        final_data_v_df[:,IndexThisVar] = ThisSeriesnew;
                    end
                end
                # print final
                if (OptionPrintFinal == "Yes") & (v == size(unique_releases,1))
                    CSV.write("Final_Dataset_For_Target.csv", final_data_v_df);
                    printstyled("Last vintage data printed in the csv file: Final_Dataset_For_Target.csv", color=:green)
                    println("")
                end
                final_data_v_array = convert(Array{Any,2},final_data_v_df[Ind_start:end,:]);
                final_forcalendar_v = convert(Array{Union{Missing,DateTime},2}, final_forcalendar_v[Ind_start:end,:]);
                # Add missings to forecasting horizon
                len_v = size(final_data_v_array[:,2:end],1);
                n_vars = size(final_data_v_array[:,2:end],2);
                data_v = missing.*ones(len_v+h, n_vars) |> Array{Union{Missing,Float64},2};
                data_v[1:len_v,:] = final_data_v_array[:,2:end];
                data_revtr[v]                 = data_v;
                date_revtr[v]                 = convert(Array{Union{Missing,DateTime},1},final_data_v_array[:,1]);
                data_revtr_year[year][v_year] = data_v;
                date_revtr_year[year][v_year] = convert(Array{Union{Missing,DateTime},1},final_data_v_array[:,1]);
                calendar_revtr[v]             = final_forcalendar_v;
                # Update vintage ID
                v += 1;
        end
    end
    return data_revtr, data_revtr_year, date_revtr, date_revtr_year, calendar_revtr;
end
