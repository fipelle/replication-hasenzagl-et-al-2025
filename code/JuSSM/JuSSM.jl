__precompile__()

module JuSSM
# ----------------------------------------------------------------------------------------------------------------------
# This module define the functions used to estimate linear state-space models with an Adaptive Metropolis-Within-Gibbs.
#
# This algorithm is structured in two blocks:
# a) The first one is a Single Component Metropolis a lá Haario et al (2003). This is used to estimate the parameters.
#    The major difference with the SCAM approach is that here the parameters are transformed so that their support is
#	 unconstrained, in order to improve mixing. Prior densities are corrected with the log-J of the transformations.
#
# b) The latent states are drawn conditional to the model parameters.
#
# ----------------------------------------------------------------------------------------------------------------------

	using DataFrames;
    using Distributions;
    using LinearAlgebra;
	using Distributed;
	using Random;

	local_path = dirname(@__FILE__);

	# -----------------------------------------------------------------------------------------------------------------
	# Types
	# -----------------------------------------------------------------------------------------------------------------

	mutable struct ParSsm{X <: Union{Missing,Float64}}
		y::Array{Union{Missing,Float64}, 2}
		d::Union{Array{X, 1}}
		Z::Union{Array{X, 1}, Array{X, 2}}
		R::Union{Array{X, 1}, Array{X, 2}}
		c::Union{Array{X, 1}}
		T::Union{Array{X, 1}, Array{X, 2}}
		Q::Union{Array{X, 1}, Array{X, 2}}
		α¹::Union{Array{X, 1}}
		P¹::Union{Array{X, 1}, Array{X, 2}}
		P̄¹::Union{Array{X, 1}, Array{X, 2}}
		λ::Array{X, 1}
		ρ::Array{X, 1}
		companion_cycle::Array{Int64, 1}
		restr_d_SPF::Bool
		logprior::X
		loglik::X
		logposterior::X
	end

	mutable struct SizeParSsm{X <: Int64}
		d::X
		Z::X
		R::X
		c::X
		T::X
		Q::X
		λ::X
		ρ::X
		θ::X
	end

	mutable struct BoolParSsm{X1 <: BitArray{1}, X2 <: BitArray{2}}
		d::X1
		Z::X2
		R::X2
		c::X1
		T::X2
		Q::X2
		λ::X1
		ρ::X1
		quarterly::X1
	end

	mutable struct PriorOpt{X <: Float64} # the logpdf for λ and ρ are constants
		N::Distributions.Normal{X}
		IG::Distributions.InverseGamma{X} #TruncatedNormal{X}
		λ::X
		ρ::X
	end

	struct output_position{X <: Int64}
		BC::X
		EP::X
		T_INFL::X
		last_GDP_state::X
		GDP::X
		INFL::X
	end

	import Base.copy
	Base.copy(m::ParSsm) 	= ParSsm([ copy(getfield(m, k)) for k = 1:length(fieldnames(typeof(m))) ]...);
    Base.copy(m::SizeParSsm) = SizeParSsm([ copy(getfield(m, k)) for k = 1:length(fieldnames(typeof(m))) ]...);
    Base.copy(m::BoolParSsm) = BoolParSsm([ copy(getfield(m, k)) for k = 1:length(fieldnames(typeof(m))) ]...);
    Base.copy(m::PriorOpt) 	= PriorOpt([ copy(getfield(m, k)) for k = 1:length(fieldnames(typeof(m))) ]...);

	# -----------------------------------------------------------------------------------------------------------------
	# Subroutines
	# -----------------------------------------------------------------------------------------------------------------

	# Estimation
	include("$local_path/subroutines/estimation/kalman_diffuse!.jl");
	include("$local_path/subroutines/estimation/JuSSM_main.jl");
	include("$local_path/subroutines/estimation/JuSSM_run.jl");

	# Julia's subroutines
	include("$local_path/subroutines/extra/ex_blkdiag.jl");
	include("$local_path/subroutines/extra/ex_inv.jl");
	include("$local_path/subroutines/extra/ex_ismember.jl");

	# Get
	include("$local_path/subroutines/get/get_logjacobian.jl");
	include("$local_path/subroutines/get/get_par_bound.jl");
	include("$local_path/subroutines/get/get_par_unb.jl");
	include("$local_path/subroutines/get/get_progress.jl");
	include("$local_path/subroutines/get/get_random_disturbance.jl");

	# Set
	include("$local_path/subroutines/set/set_par.jl");


	# -----------------------------------------------------------------------------------------------------------------
	# Export functions
	# -----------------------------------------------------------------------------------------------------------------

	export kalman_diffuse!, JuSSM_main, ex_blkdiag, ex_inv, ex_ismember, set_par_fast!, ParSsm, SizeParSsm, BoolParSsm, PriorOpt, output_position;
end
