function ex_inv(x)
# ----------------------------------------------------------------------------------------------------------------------
# Custom inv function
# ----------------------------------------------------------------------------------------------------------------------

    xInv = x\Matrix{Float64}(I, size(x)[1], size(x)[2]);
    #eye(size(x)[1], size(x)[2]);

    return xInv;
end
