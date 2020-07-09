function get_par_unb(θ_bound::X, MIN::X, MAX::X, opt_transf::Array{Int64, 1}) where X<:Array{Float64,1}
# ----------------------------------------------------------------------------------------------------------------------
# Change domain: from bounded to unbounded domain
# - Transformation: 1 natural logarithm, 2 no transformations, 3 generalized logit
#
# Author: Filippo Pellegrino, f.pellegrino1@lse.ac.uk
# ----------------------------------------------------------------------------------------------------------------------

     i1 = opt_transf .== 1;
     i2 = opt_transf .== 2;
     i3 = opt_transf .== 3;

     θ_unb     = zeros(size(θ_bound));
     θ_unb[i1] = log.(θ_bound[i1]-MIN[i1]);
     θ_unb[i2] = θ_bound[i2];
     θ_unb[i3] = log.((θ_bound[i3]-MIN[i3])./(MAX[i3]-θ_bound[i3]));

     return θ_unb;
end
