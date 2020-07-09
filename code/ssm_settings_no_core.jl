function ssm_settings(y, h, nDraws, burnin, σʸ, quarterly_position, estim_bool, ind_restr_states; size_only_bool=false, oos_details="")
# ----------------------------------------------------------------------------------------------------------------------
# Define the basic structure of the state-space parameters
# ----------------------------------------------------------------------------------------------------------------------

    n  = size(y)[2];
    nR = length(ind_restr_states);


    # -----------------------------------------------------------------------------------------------------------------
    # Observation equations
    # -----------------------------------------------------------------------------------------------------------------

    d    = zeros(n+nR);
    d[2] = 1;
    Z    = [[ones(4); 0; ones(n-5)] [zeros(1,3); ones(3,3); zeros(1,3); ones(n-5,3)] [zeros(4); ones(n-4)] zeros(n) [zeros(5); 1 ./ σʸ[end-2:end]]];

    # Other loadings
    Z1 = permutedims(kron(Matrix(I,3,3), [1, 0, 1]));                                                                                     # idio C, idio C+, idio trend
    Z2 = permutedims(kron(Matrix(I,2,2), [1, 0]));                                                                                        # idio C, idio C+
    Z  = [Z cat(dims=[1,2], [1 0 1 ./ σʸ[1]], [1.0 0.0], Z1, Z2, [1.0 0.0 1.0])]; # No idio trend for SPF

    # Link GDP SPF with GDP (trend and drift)
    Z[2, 10] = 3 ./ σʸ[2];

    # Mixed frequency restriction for quarterly variables
    nQ                                        = convert(Int64, sum(quarterly_position));    # No. quarterly variables
    Z                                         = [Z zeros(n, 2*nQ)];                         # C(-1), C(-2), T(-1), T(-2)
    Z[quarterly_position.==1, end-2*nQ+1:end] = kron(Matrix(I,nQ,nQ), ones(2, 1))';                 # Link C(-1), C(-2), T(-1), T(-2)

    # Set linear restrictions in Z
    for i=1:nR

        # Index
        ind_var_restr = findall(Z[:, ind_restr_states[i]])[1];
        first_not_na  = findall(ismissing.(y)[:, ind_var_restr] .== false)[1];

        # Apply restriction
        ZRᵢ = zeros(1, size(Z)[2]);
        ZRᵢ[1, ind_restr_states[i]] = 1;
        Z = [Z; ZRᵢ];

        # Generate auxiliary variables for the linear restrictions
        if nR > 0

            # Generate
            auxiliary = missing.*ones(size(y)[1]) |> Array{Union{Float64, Missing}};
            auxiliary[first_not_na] = 0.0;

            # Append
            y = [y auxiliary];
        end

        @info("Variable $ind_var_restr (state $(ind_restr_states[i])) restricted at t$first_not_na");
    end

    # Irregular components
    R = cat(dims=[1,2], 1e-4.*Matrix(I,n,n), zeros(nR, nR));

    # Indeces for observation equations
    d_ind = d .!= 0;
    Z_ind = zeros(size(Z)) .!= 0;
    R_ind = zeros(size(R)) .!= 0; # R .!= 0;

    # Projections
    Z_ind[2:4, [1,2,3,4]] .= true; # Real               -> PC cycle, t, t-1, t-2, t-3
    Z_ind[6:n, [1,2,3,4]] .= true; # Prices             -> PC cycle, t, t-1, t-2, t-3
    Z_ind[6:n, 5]         .= true; # Prices             -> EP cycle, t (no lags)

    # estim_bool == true is for consumer expectations in percentage of balance (not the default option)
    if estim_bool == true
        Z_ind[n, 7] .= true; # EC expectations -> T
    end

    restr_d_SPF = false;


    # -----------------------------------------------------------------------------------------------------------------
    # Transition equations
    # -----------------------------------------------------------------------------------------------------------------

    c             = zeros(size(Z)[2]);
    ind_trends    = [10; 18]; # GDP DRIFT, EMPL
    c[ind_trends] .= 1;  # random walk drift

    T_c     = convert(Array{Float64, 2}, [1 0; 0 0]);
    T_ct    = convert(Array{Float64, 2}, [1 0 0; 0 0 0; 0 0 1]);
    T_c_ext = convert(Array{Float64, 2}, [1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 1 0]);
    Q_c_ext = convert(Array{Float64, 2}, [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);

    T = cat(dims=[1,2], T_c_ext, [T_ct for i=1:2]..., T_c, [T_ct for i=1:3]..., [T_c for i=1:2]..., T_ct, zeros(2*nQ, 2*nQ));
    Q = cat(dims=[1,2], Q_c_ext, [T_ct for i=1:2]..., T_c, [T_ct for i=1:3]..., [T_c for i=1:2]..., T_ct, zeros(2*nQ, 2*nQ));

    # Indeces for transition equations
    c_ind = c .!= 0;
    T_ind = zeros(size(T)) .== 1;
    Q_ind = Q .== 1;

    # Initial conditions for the non-stationary states
    P̄_c   = convert(Array{Float64, 2}, [0 0; 0 0]);
    P̄_ct  = convert(Array{Float64, 2}, [0 0 0; 0 0 0; 0 0 1]);
    P̄¹    = cat(dims=[1,2], zeros(4,4), [P̄_ct for i=1:2]..., P̄_c, [P̄_ct for i=1:3]..., [P̄_c for i=1:2]..., P̄_ct, zeros(2*nQ, 2*nQ));

    # Initial conditions
    α¹ = zeros(size(c));
    P¹ = zeros(size(P̄¹));

    # Trigonometric states
    λ_c   = convert(Array{Float64, 1}, [1; 0]);
    λ_ct  = convert(Array{Float64, 1}, [1; 0; 0]);
    λ     = vcat([1; zeros(3)], [λ_ct for i=1:2]..., λ_c, [λ_ct for i=1:3]..., [λ_c for i=1:2]..., λ_ct, zeros(2*nQ));
    ρ     = copy(λ);
    λ_ind = λ .!= 0;
    ρ_ind = copy(λ_ind);

    companion_cycle = copy(λ_ind) |> Array{Int64,1};
    companion_cycle[1] = 3;


    # -----------------------------------------------------------------------------------------------------------------
    # Mixed-frequency
    # -----------------------------------------------------------------------------------------------------------------

    ind_Q = findall(quarterly_position.==1);

    for i=1:nQ
        ind_Zᵢ = findall(Z[ind_Q[i], :] .!= 0)[1:end-2];
        ind_Tᵢ = size(T)[1]-2*nQ+1+(i-1)*2;

        T[ind_Tᵢ, ind_Zᵢ]  .= 1.0; # lagged forecast
        T[ind_Tᵢ+1, ind_Tᵢ] = 1.0; # lagged state
    end


    # -----------------------------------------------------------------------------------------------------------------
    # JuSSM
    # -----------------------------------------------------------------------------------------------------------------

    par_ind = BoolParSsm(d_ind, Z_ind, R_ind, c_ind, T_ind, Q_ind, λ_ind, ρ_ind, quarterly_position.==1);
    par     = ParSsm(permutedims(y), d, Z, R, c, T, Q, α¹, P¹, P̄¹, λ, ρ, companion_cycle, restr_d_SPF, 0.0, 0.0, 0.0);

    if size_only_bool == false

        distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, par, par_size, distr_par =
            JuSSM_main(par, h, nDraws, burnin, par_ind, oos_details=oos_details);

        return distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, par, par_ind, par_size, distr_par;

    else
        size_θ = sum(par_ind.d) + sum(sum(par_ind.Z)) + sum(sum(par_ind.R)) +
                 sum(par_ind.c) + sum(sum(par_ind.T)) + sum(sum(par_ind.Q)) +
                 sum(sum(par_ind.λ)) + sum(sum(par_ind.ρ));

        return size(Z)[2], size_θ;
    end
end
