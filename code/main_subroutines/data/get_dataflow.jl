function get_dataflow(data, date, data_path::String, h::Int64, oos_start_date::DateTime)

    #=
    --------------------------------------------------------------------------------------------------------------------
    Description:   Get data vintages from the excel file. Dates must be in "date" format in data file.
    Output:        a) Data vintages, b) Data vintages, per year, c) years, d) releases_per_year
    Author:        Filippo Pellegrino
    Email:         filippopellegrino@gmail.com
    --------------------------------------------------------------------------------------------------------------------
    =#

    #=
    -----------------------------------
    Initialise
    -----------------------------------
    =#

    # Size
    n = size(data)[2];

    # Load data from the excel file
    raw_calendar = DataFrame(XLSX.readtable(data_path, "calendar")...);
    reference    = vcat([names(raw_calendar)[1+i] |> String |> DateTime for i=1:size(raw_calendar,2)-1]...);
    releases     = raw_calendar[:, 2:end] |> Array{Union{Missing,DateTime},2};

    # Unique releases and years
    unique_releases = sort(unique(skipmissing(releases[:])));
    years           = vcat([Dates.year(unique_releases[i]) for i=1:length(unique_releases)]...);
    unique_years    = sort(unique(years));

    # Number of releases per year
    releases_per_year = vcat([sum(years .== unique_years[i]) for i=1:length(unique_years)]...);

    # Arrays of data vintages
    data_vintages      = Array{Array{Union{Missing,Float64},2}}(undef,length(unique_releases),1);
    data_vintages_year = Array{
        Array{
            Array{Union{Missing,Float64},2}
            ,1}}(undef,length(unique_years));
    for i=1:length(unique_years)
        data_vintages_year[i] = Array{Array{Union{Missing,Float64},2},1}(undef,releases_per_year[i]);
    end

    # Error management: remove first "discard" vintages from data_vintages
    check_discard = true;
    len_0         = 1;
    discard       = 0;

    # NA -> NaN
    releases[ismissing.(releases)] .= DateTime("01-01-1500", "dd-mm-yyyy");


    #=
    -----------------------------------
    Build collections of data vintages
    -----------------------------------
    =#

    # Vintage ID
    v = 1;

    for year=1:length(unique_years)
        for v_year=1:releases_per_year[year]

            #@info(year, v_year, v)
            #=
            -----------------------------------
            Indeces
            -----------------------------------
            =#

            # Find the release in "releases"
            bool_release           = releases .== unique_releases[v];
            releases_per_ref       = sum(Array{Int64,2}(bool_release), dims=1)'[:,1];
            bool_release           = bool_release[:, findall(releases_per_ref .> 0)];

            # Reference period
            reference_v = reference[findall(releases_per_ref .> 0)];

            # Vintage v
            if v == 1
                len_v = findall(date.==maximum(reference_v))[1];
            else
                data_v_1 = data_vintages[v-1];
                len_v    = max(findall(date.==maximum(reference_v))[1], size(data_v_1)[1]-h);
            end

            data_v = missing.*ones(len_v+h, n) |> Array{Union{Float64,Missing}};

            #=
            -----------------------------------
            Data up to previous release
            -----------------------------------
            =#

            # First release
            if v == 1
                data_v[1:len_v-1, :] = copy(data[1:len_v-1, :]);
                len_0                = len_v-1;

            # Subsequent releases
            else
                data_v[1:size(data_v_1)[1], :] = copy(data_v_1);
            end

            #=
            -----------------------------------
            Output
            -----------------------------------
            =#

            for i=1:size(bool_release, 2)

                # Row and column corresponding to the new release
                data_row = findall(date.==reference_v[i])[1];
                data_col = findall(bool_release[:, i]);

                # Update data_v
                data_v[data_row, data_col] = copy(data[data_row, data_col]);
            end

            # Error management: remove first "discard" vintages from data_vintages
            if (check_discard == true) && (v > 1)
                if sum(ismissing.(data_v[len_0+1, :])) == 0
                    discard       = copy(v-1);
                    check_discard = false;
                end
            end

            data_vintages[v]                 = data_v;
            data_vintages_year[year][v_year] = data_v;

            # Update vintage ID
            v += 1;
        end
    end

    #=
    -----------------------------------
    Remove first "discard" vintages
    -----------------------------------
    =#

    if unique_releases[discard+1] > oos_start_date
        error("For each variable we need at least one release date before the OOS start date");
    elseif unique_releases[discard+1] < oos_start_date
        discard = findfirst(unique_releases.>=oos_start_date)-1;
    end

    @info("get_dataflow > OOS Start Date: $(unique_releases[discard+1])")

    # Full
    data_vintages = data_vintages[discard+1:end];

    # Per year
    years_to_remove       = cumsum(releases_per_year) .< discard + 1;
    remainder             = discard - sum(releases_per_year[years_to_remove]);
    data_vintages_year    = data_vintages_year[years_to_remove.==false];
    data_vintages_year[1] = data_vintages_year[1][remainder+1:end];
    unique_years          = unique_years[years_to_remove.==false];
    releases_per_year     = vcat([length(data_vintages_year[i]) for i=1:length(data_vintages_year)]...);

    return data_vintages, data_vintages_year, unique_years, releases_per_year;
end
