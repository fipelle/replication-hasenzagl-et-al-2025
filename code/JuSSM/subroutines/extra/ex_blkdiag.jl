function ex_blkdiag(args::Union{Array{W}, W}...) where W <: Union{Missing,Float64}

    #=
    --------------------------------------------------------------------------------------------------------------------
    Description:   Custom built-in function: blkdiag
    --------------------------------------------------------------------------------------------------------------------
    =#

    out = args[1];

    for arg in args[2:end]
        out = blkdiag_fun(out, arg);
    end

    return out
end


function blkdiag_fun(x::Union{Array{W}, W}, y::Union{Array{W}, W}) where W <: Union{Missing,Float64};

    #=
    --------------------------------------------------------------------------------------------------------------------
    Description:   Custom built-in function: blkdiag, single iteration
    --------------------------------------------------------------------------------------------------------------------
    =#

    # Initialise
    if ndims(x) > 1
        rx, cx = size(x);
    else
        rx = size(x)[1];
        cx = 1;
    end

    if ndims(y) > 1
        ry, cy = size(y);
    else
        ry = size(y)[1];
        cy = 1;
    end

    # Blkdiag
    out                         = zeros(rx+ry, cx+cy);
    out[1:rx, 1:cx]             = x;
    out[rx+1:rx+ry, cx+1:cx+cy] = y;

    return out;
end
