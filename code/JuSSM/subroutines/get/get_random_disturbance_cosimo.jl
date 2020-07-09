function get_random_disturbance(x, Σ)

     # Init random perturbation
     rnd_per     = zeros(x);
     rnd_per_ind = (sum(Σ, dims = 2).!=0)[:,1];

     if sum(rnd_per_ind) > 0

          # Get random perturbation
          rnd_per_cov                = Σ[rnd_per_ind, rnd_per_ind];
          rnd_per[rnd_per_ind] = rand(MvNormal(zeros(sum(rnd_per_ind)), rnd_per_cov));
     end

     return rnd_per;
end
