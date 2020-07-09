function get_par_bound(θ_unb::X, MIN::X, MAX::X, opt_transf::Array{Int64, 1}) where X<:Array{Float64,1}
# ----------------------------------------------------------------------------------------------------------------------
# Change domain: from bounded to unbounded domain
# - Transformation: 1 natural logarithm, 2 no transformations, 3 generalized logit
#
# Author: Filippo Pellegrino, f.pellegrino1@lse.ac.uk
# ----------------------------------------------------------------------------------------------------------------------

     i1 = opt_transf .== 1;
     i2 = opt_transf .== 2;
     i3 = opt_transf .== 3;

     θ_bound     = zeros(size(θ_unb));
     θ_bound[i1] = exp.(θ_unb[i1])+MIN[i1];
     θ_bound[i2] = θ_unb[i2];
     θ_bound[i3] = (MIN[i3] .+ MAX[i3].*exp.(θ_unb[i3]))./(1 .+ exp.(θ_unb[i3]));

     return θ_bound;
end
