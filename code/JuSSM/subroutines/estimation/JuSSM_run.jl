function JuSSM_run(θ_unb::Array{Float64,1}, par::ParSsm, h::Int64, par_ind::BoolParSsm, par_size::SizeParSsm,
                   prior_opt::PriorOpt, MIN::Array{Float64, 1}, MAX::Array{Float64, 1}, opt_transf::Array{Int64, 1},
                   nDraws::Int64, burnin::Int64, algorithm_name::String, oos_details::String, gibbs_step=1::Int64)

    # ------------------------------------------------------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------------------------------------------------------

    # Pre-allocate memory
    k             = size(par.T)[1];
    n, m          = size(par.y);
    distr_α       = zeros(k, m, nDraws-burnin);
    distr_fcst    = zeros(m, n, nDraws-burnin);
    distr_par     = Array{Any}(undef, nDraws-burnin);
    chain_θ_unb   = zeros(par_size.θ, nDraws);
    chain_θ_bound = zeros(par_size.θ, nDraws);

    # Initial values
    θ_bound           = get_par_bound(θ_unb, MIN, MAX, opt_transf);
    par.logposterior  = -Inf;
    par_prop          = copy(par);
    apriori_rejection = [0.0];

    # Mwg scale
    mwg_scale = ones(par_size.θ);

    if oos_details != ""
        @info("$algorithm_name ($oos_details) > Initialisation, $(par_size.θ) parameters to estimate");
    else
        @info("$algorithm_name > Initialisation, $(par_size.θ) parameters to estimate");
    end


    # ------------------------------------------------------------------------------------------------------------------
    # Start estimation
    # ------------------------------------------------------------------------------------------------------------------

    for draw=1:nDraws

        last_draw = draw-1;

        if (draw == 11) || ((draw > 1) && mod(last_draw, 100) == 0)

            last_logposterior = round(par.logposterior, digits = 2);

            if oos_details != ""
                print("$algorithm_name ($oos_details) > Progress: $last_draw/$nDraws draws | Logposterior: $last_logposterior");
            else
                print("$algorithm_name > Progress: $last_draw/$nDraws draws | Logposterior: $last_logposterior");
            end

            get_progress(chain_θ_unb, 1, last_draw, 1);

            if last_draw > burnin+100
                get_progress(chain_θ_unb, burnin+1, last_draw, 1, "Acceptance rate (PB)");
            end

            @info("");
        end


        # --------------------------------------------------------------------------------------------------------------
        # Block 1: Propose θ* | θ_{draw-1}
        # --------------------------------------------------------------------------------------------------------------

        # Component-wise algorithm (shuffle the order of the parameters for every draw)
        for par_id = shuffle(1:par_size.θ)

            # Variant of SCAM (see working paper)
            # Note: optimal acceptance rate for single component algorithms is 44%
            if last_draw > 10
                acc_rate           = get_progress(chain_θ_unb[par_id, :], 1, last_draw, 0);
                mwg_scale[par_id] *= exp.((acc_rate-44.0)/100);
            end

            # Draw new candidate
            θ_prop_unb          = copy(θ_unb);
            θ_prop_unb[par_id] += mwg_scale[par_id]*randn()
            θ_prop_bound        = get_par_bound(θ_prop_unb, MIN, MAX, opt_transf);

            set_par!(θ_prop_bound, θ_prop_unb, par_prop, opt_transf, MIN, MAX, par_ind, par_size, prior_opt, apriori_rejection);

            if apriori_rejection[1] == 0

                # Accept / Reject
                if rand(Uniform(0, 1)) < exp.(par_prop.logposterior-par.logposterior)

                    # Accept the proposed state-space parameters
                    par     = copy(par_prop);
                    θ_unb   = copy(θ_prop_unb);
                    θ_bound = copy(θ_prop_bound);
                end
            end
        end

        # Chains
        chain_θ_unb[:, draw]   = θ_unb;
        chain_θ_bound[:, draw] = θ_bound;


        # --------------------------------------------------------------------------------------------------------------
        # Block 2: Draw alpha (Gibbs step)
        # --------------------------------------------------------------------------------------------------------------

        if draw > burnin && gibbs_step == 1

            # Draw α conditional to the state-space parameters
            α_draw, _                  = kalman_diffuse!(par, 0, 1, 1);
            distr_α[:, :, draw-burnin] = α_draw;
            distr_par[draw-burnin]     = par;

            # Forecast
            distr_fcst[:, :, draw-burnin] = (par.Z*α_draw)';
        end
    end

    return chain_θ_unb, chain_θ_bound, distr_α, distr_fcst, par, distr_par;
end
