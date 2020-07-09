function set_par!(θ_bound, θ_unb, par, opt_transf, MIN, MAX, par_ind, par_size, prior_opt, apriori_rejection)

# ----------------------------------------------------------------------------------------------------------------------
# Update the parameters of the state-space model
#
# Author: Filippo Pellegrino, f.pellegrino1@lse.ac.uk
# ----------------------------------------------------------------------------------------------------------------------


     if sum(ismissing.(θ_bound)) > 0 || sum(isinf.(θ_bound)) > 0 || isreal(θ_bound) != 1
          apriori_rejection[1] = 1;

     else

          # ------------------------------------------------------------------------------------------------------------
          # Set state-space parameters and par.logprior
          # ------------------------------------------------------------------------------------------------------------

          iend         = 0;
          par.logprior = 0;

          # Observation equations

          if par_size.R > 0
               par.R[par_ind.R .== true] = θ_bound[1:par_size.R];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.IG, par.R[par_ind.R .== true]));
               iend                      = par_size.R;
          end

          if par_size.d > 0
               par.d[par_ind.d .== true] = θ_bound[iend+1:iend+par_size.d];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.N, par.d[par_ind.d .== true]));
               iend                      = iend+par_size.d;
          end

          if par_size.Z > 0
               par.Z[par_ind.Z .== true] = θ_bound[iend+1:iend+par_size.Z];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.N, par.Z[par_ind.Z .== true]));
               iend                      = iend+par_size.Z;
          end


          # Transition equations

          if par_size.Q > 0
               par.Q[par_ind.Q .== true] = θ_bound[iend+1:iend+par_size.Q];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.IG, par.Q[par_ind.Q .== true]));
               iend                      = iend+par_size.Q;
          end

          if par_size.c > 0
               par.c[par_ind.c .== true] = θ_bound[iend+1:iend+par_size.c];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.N, par.c[par_ind.c .== true]));

               if par.restr_d_SPF == true
                   par.d[2] = θ_bound[iend+1]*(12+11+10);
               end

               iend = iend+par_size.c;
          end

          if par_size.T > 0 || par_size.λ > 0 || par_size.ρ > 0

               # Set T
               par.T[par_ind.T .== true] = θ_bound[iend+1:iend+par_size.T];
               par.logprior              = par.logprior + sum(logpdf.(prior_opt.N, par.T[par_ind.T .== true]));
               iend                      = iend+par_size.T;

               # Trigonometric states: update T, λ and ρ. Adjust Q and P¹
               if par_size.λ > 0 || par_size.ρ > 0
                    par.λ[par_ind.λ .== true] = θ_bound[iend+1:iend+par_size.λ];
                    par.ρ[par_ind.ρ .== true] = θ_bound[iend+par_size.λ+1:end];
                    par.logprior              = par.logprior + prior_opt.λ + prior_opt.ρ;

                    bool_trig = ((par_ind.λ .== true) .+ (par_ind.ρ .== true)) .>= 1;
                    find_trig = findall(bool_trig);

                    for i=1:sum(bool_trig)
                         j = find_trig[i];

                         par.T[j:j+1, j:j+1]  = par.ρ[j]*[[cos(par.λ[j]) sin(par.λ[j])]; [-sin(par.λ[j]) cos(par.λ[j])]];
                         par.Q[j+1, j+1]      = par.Q[j, j];
                         par.P¹[j:j+par.companion_cycle[j], j:j+par.companion_cycle[j]] .= (par.Q[j, j] ./ (1-par.ρ[j]^2)) .* Matrix(I, par.companion_cycle[j]+1, par.companion_cycle[j]+1);
                    end
               end
          end


          # ------------------------------------------------------------------------------------------------------------
          # Mixed frequency
          # ------------------------------------------------------------------------------------------------------------

          nQ    = sum(par_ind.quarterly);
          ind_Q = findall(par_ind.quarterly);

          for i=1:nQ
               ind_Zᵢ = findall(par.Z[ind_Q[i], :] .!= 0)[1:end-2];
               ind_Tᵢ = size(par.T)[1]-2*nQ+1+(i-1)*2;

               par.T[ind_Tᵢ, ind_Zᵢ] = copy(par.Z[ind_Q[i], ind_Zᵢ]); # lagged forecast
          end


          # ------------------------------------------------------------------------------------------------------------
          # par.loglikelihood and par.logposterior
          # ------------------------------------------------------------------------------------------------------------

          kalman_diffuse!(par, 1, 0, 0); # estimate loglikelihood
          par.logprior     = par.logprior + sum(get_logjacobian(θ_unb, MIN, MAX, opt_transf));
          par.logposterior = par.loglik + par.logprior;
     end
end
