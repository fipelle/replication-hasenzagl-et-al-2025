include("./monthly2quarterly.jl");

function read_data(data_path)

     date_monthly   = [];
     data_monthly   = [];
     date_quarterly = [];
     data_quarterly = [];

     # Monthly data
     try
          f_monthly    = DataFrame(XLSX.readtable(data_path, "monthly")...);
          date_monthly = f_monthly[18:end, 2] |> Array{Date,1};
          data_monthly = f_monthly[18:end, 3:end] |> Array{Union{Float64, Missing},2};

     catch
          print("There are no monthly variables in $data_path \n");
     end

     # Quarterly data
     try
          f_quarterly    = DataFrame(XLSX.readtable(data_path, "quarterly")...);
          date_quarterly = f_quarterly[18:end, 2] |> Array{Date,1};;
          data_quarterly = f_quarterly[18:end, 3:end] |> Array{Union{Float64, Missing},2};

     catch
          print("There are no quarterly variables in $data_path \n");
     end

     # Merge
     if size(data_monthly)[1] == 0 && size(data_quarterly)[1] > 0
          date = date_quarterly;
          data = data_quarterly;
          nM   = 0;
          nQ   = size(data)[2];

     elseif size(data_monthly)[1] > 0 && size(data_quarterly)[1] == 0
          date = date_monthly;
          data = data_monthly;
          nM   = size(data)[2];
          nQ   = 0;

     elseif size(data_monthly)[1] == 0 && size(data_quarterly)[1] == 0
          date     = [];
          data     = [];
          nM       = 0;
          nQ       = 0;
          MNEMONIC = [];
          print("There is no data in $data_path \n");

     else
          Xq   = monthly2quarterly(data_quarterly, size(data_monthly)[1]);
          date = date_monthly;
          data = [[data_monthly; missing.*ones(size(Xq,1)-size(data_monthly,1), size(data_monthly,2))] Xq];
          nM   = size(data_monthly)[2];
          nQ   = size(data_quarterly)[2];
     end

     data = data |> Array{Union{Float64, Missing},2};

     # MNEMONIC
     if nM != 0 || nQ != 0
         try
              info_data = DataFrame(XLSX.readtable(data_path, "transf")...);
              MNEMONIC  = info_data[1:end, 2] |> Array{String,1};

          catch
              print("There are no mnemonics in $data_path \n");
              MNEMONIC = [];
          end
     end

     return data, date, nM, nQ, MNEMONIC;
end

function read_data_info(data_path)

    # Load Excel file
    f_info = DataFrame(XLSX.readtable(data_path, "info")...);

    # Mnemonic
    MNEMONIC = f_info[!,:MNEMONIC] |> Array{String,1};

    # nM and nQ
    frequency = f_info[!,:FREQ] |> Array{String,1};
    nM = sum(frequency .== "m");
    nQ = sum(frequency .== "q");

    # transf and transf_annex
    transf = f_info[!,:TRANSF] |> Array{Int64,1};
    transf_arg1 = f_info[!,:TRANSF_ARG1] |> Array{Any,1};
    transf_arg2 = f_info[!,:TRANSF_ARG2] |> Array{Any,1};

    # Return output
    return MNEMONIC, nM, nQ, transf, transf_arg1, transf_arg2;
end
