function ex_inv(x)
# ----------------------------------------------------------------------------------------------------------------------
# Custom inv function
#
# Author: Filippo Pellegrino, f.pellegrino1@lse.ac.uk
# ----------------------------------------------------------------------------------------------------------------------

    xInv = x\Matrix{Float64}(I, size(x)[1], size(x)[2]);
    #eye(size(x)[1], size(x)[2]);

    return xInv;
end
