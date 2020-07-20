function JuSSM_main(par::ParSsm, h::Int64, nDraws::Int64, burnin::Int64, par_ind::BoolParSsm; oos_details::String="")
# ----------------------------------------------------------------------------------------------------------------------
# JuSSM: mainframe
#
# ----------------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    # Initialise
    # -----------------------------------------------------------------------------------------------------------------

    # Dimensions
    k    = size(par.T)[1];
    n, m = size(par.y);

    # No. parameters to be estimated
    par_size = SizeParSsm(sum(par_ind.d),
                                sum(sum(par_ind.Z)),
                                sum(sum(par_ind.R)),
                                sum(par_ind.c),
                                sum(sum(par_ind.T)),
                                sum(sum(par_ind.Q)),
                                sum(sum(par_ind.λ)),
                                sum(sum(par_ind.ρ)),
                                sum(par_ind.d) + sum(sum(par_ind.Z)) + sum(sum(par_ind.R)) +
                                    sum(par_ind.c) + sum(sum(par_ind.T)) + sum(sum(par_ind.Q)) +
                                    sum(sum(par_ind.λ)) + sum(sum(par_ind.ρ)));


    # -----------------------------------------------------------------------------------------------------------------
    # Lower and upper bounds
    # -----------------------------------------------------------------------------------------------------------------

    # Numerical constants
    xi        = 1e-3; # it must be small
    MIN_var   = 0;
    MIN_coeff = -Inf;
    MIN_λ     = xi;
    MIN_ρ     = xi;
    MAX_var   = Inf;
    MAX_coeff = Inf;
    MAX_λ     = pi;
    MAX_ρ     = 0.98;

    # MIN
    MIN = [MIN_var*ones(par_size.R);
                 MIN_coeff*ones(par_size.d + par_size.Z);
                 MIN_var*ones(par_size.Q);
                 MIN_coeff*ones(par_size.c + par_size.T);
                 MIN_λ*ones(par_size.λ);
                 MIN_ρ*ones(par_size.ρ)];

    # MAX
    MAX = [MAX_var*ones(par_size.R);
                 MAX_coeff*ones(par_size.d + par_size.Z);
                 MAX_var*ones(par_size.Q);
                 MAX_coeff*ones(par_size.c + par_size.T);
                 MAX_λ*ones(par_size.λ);
                 MAX_ρ*ones(par_size.ρ)];


    # -----------------------------------------------------------------------------------------------------------------
    # Set prior parameters
    # -----------------------------------------------------------------------------------------------------------------

    # Define prior distribution objects
    prior_opt = PriorOpt(Normal(0, 1/xi),
                                InverseGamma(3, 1), #TruncatedNormal(0, 1, 0, Inf),
                                par_size.λ*logpdf.(Uniform(MIN_λ, MAX_λ), MIN_λ),
                                par_size.ρ*logpdf.(Uniform(MIN_ρ, MAX_ρ), MIN_ρ));

    # Transformations: 1 natural logarithm, 2 no transformations, 3 generalized logit
    opt_transf = convert(Array{Int64, 1}, [1*ones(par_size.R);
                                                 2*ones(par_size.d + par_size.Z);
                                                 1*ones(par_size.Q);
                                                 2*ones(par_size.c + par_size.T);
                                                 3*ones(par_size.λ);
                                                 3*ones(par_size.ρ)]);


    # -----------------------------------------------------------------------------------------------------------------
    # Initialise JuSSM_run
    # -----------------------------------------------------------------------------------------------------------------

    # Initial guess for θ_ini_bound
    θ_ini_bound = [ones(par_size.R);
                         zeros(par_size.d);
                         ones(par_size.Z);
                         ones(par_size.Q);
                         zeros(par_size.c);
                         ones(par_size.T);
                         (2*pi/96)*ones(par_size.λ);
                         0.9*ones(par_size.ρ)]; # As in Harvey and Trimbur 2003

    # Run init
    θ_ini_unb = get_par_unb(θ_ini_bound, MIN, MAX, opt_transf);


    # -----------------------------------------------------------------------------------------------------------------
    # JuSSM_run
    # -----------------------------------------------------------------------------------------------------------------

    algorithm_name = "JuSSM";
    chain_θ_unb    = Array{Union{Missing,Float64}}(undef, 1);
    chain_θ_bound  = Array{Union{Missing,Float64}}(undef, 1);
    distr_α        = Array{Union{Missing,Float64}}(undef, 1);
    distr_fcst     = Array{Union{Missing,Float64}}(undef, 1);
    distr_par      = Array{Union{Missing,Float64}}(undef, 1);

    chain_θ_unb, chain_θ_bound, distr_α, distr_fcst, par, distr_par =
        JuSSM_run(θ_ini_unb, par, h, par_ind, par_size, prior_opt, MIN, MAX, opt_transf, nDraws, burnin, algorithm_name,
                oos_details, 1);

    return distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, par, par_size, distr_par;
end
