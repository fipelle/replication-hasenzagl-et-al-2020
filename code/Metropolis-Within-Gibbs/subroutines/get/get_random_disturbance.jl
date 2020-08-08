#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function get_random_disturbance(x, Σ)

     # Init random perturbation
     rnd_per     = zeros(x);
     rnd_per_ind = sum(Σ, dims=2)[:] .!= 0;

     if sum(rnd_per_ind) > 0

          # Get random perturbation
          rnd_per_cov                   = Σ[findall(rnd_per_ind), findall(rnd_per_ind)];
          rnd_per[findall(rnd_per_ind)] = rand(MvNormal(zeros(sum(rnd_per_ind)), rnd_per_cov));
     end

     return rnd_per;
end
