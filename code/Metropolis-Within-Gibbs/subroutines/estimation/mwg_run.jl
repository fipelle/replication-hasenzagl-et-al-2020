#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function mwg_run(θ_unb::Array{Float64,1}, par::ParSsm, h::Int64, par_ind::BoolParSsm, par_size::SizeParSsm,
                  prior_opt::PriorOpt, MIN::Array{Float64, 1}, MAX::Array{Float64, 1}, opt_transf::Array{Int64, 1},
                  nDraws::Int64, burnin::Int64, mwg_const::Float64, algorithm_name::String, t::Int64, end_oos::Int64,
                  Σ::Array{Float64, 2}, gibbs_step::Int64)

     # -----------------------------------------------------------------------------------------------------------------
     # Initialisation
     # -----------------------------------------------------------------------------------------------------------------

     # Pre-allocate memory
     k             = size(par.T)[1];
     n, m          = size(par.y);
     distr_α       = zeros(k, m, nDraws-burnin);
     distr_fcst    = zeros(h, n, nDraws-burnin);
     distr_par     = Array{Any}(undef, nDraws-burnin);
     chain_θ_unb   = zeros(par_size.θ, nDraws);
     chain_θ_bound = zeros(par_size.θ, nDraws);

     # Initial values
     θ_bound           = get_par_bound(θ_unb, MIN, MAX, opt_transf);
     par.logposterior  = -Inf;
     par_prop          = copy(par);
     apriori_rejection = [0.0];

     # Parallel computing: random jumps for the proposal distribution
     n_workers        = length(workers());
     nDraws_worker    = convert(Int64, floor(nDraws/n_workers));
     θ_unb_jumps_pmap = pmap(get_mwg_jump,
                             [Σ for i=1:n_workers],
                             [mwg_const for i=1:n_workers],
                             [nDraws_worker for i=1:n_workers]);

     θ_unb_jumps = reduce(hcat, θ_unb_jumps_pmap);

     # Additional constants
     workers_vec = repeat(collect(1:n_workers), convert(Int64, floor(100/n_workers)), 1);
     n_fields    = length(fieldnames(typeof(par)));
     n_iter      = convert(Int64, ceil(n_fields/n_workers));


     # -----------------------------------------------------------------------------------------------------------------
     # Start estimation
     # -----------------------------------------------------------------------------------------------------------------

     for draw=1:nDraws

          last_draw = draw-1;

          if draw > 1 && mod(last_draw, 1000) == 0

               last_logposterior = round(par.logposterior, digits=2);

               # Display the estimation status
               if t != 0 && end_oos != 0
                    print("OOS > Running the $t-th iteration (out of $end_oos) \n");
               end

               print("$algorithm_name > Executed $last_draw draws (out of $nDraws): \n");
               print("- Logposterior: $last_logposterior\n")
               get_progress(chain_θ_unb, 1, last_draw, 1);
               if last_draw > burnin+100
                    get_progress(chain_θ_unb, burnin+1, last_draw, 1, "Acceptance rate (post burn-in)");
               end
               print("\n");
          end


          # ------------------------------------------------------------------------------------------------------------
          # Block 1: Propose θ* | θ_{draw-1}
          # ------------------------------------------------------------------------------------------------------------

          θ_prop_unb   = θ_unb + θ_unb_jumps[:, draw];
          θ_prop_bound = get_par_bound(θ_prop_unb, MIN, MAX, opt_transf);

          set_par!(θ_prop_bound, θ_prop_unb, par_prop, opt_transf, MIN, MAX, par_ind, par_size, prior_opt, apriori_rejection);

          if apriori_rejection[1] == 0

               # Accept / Reject
               if rand(Uniform(0, 1)) < exp(par_prop.logposterior-par.logposterior)

                    # Accept the proposed state-space parameters
                    # - Copy the fields of par_prop in par (in parallel)
                    begin
                         par_list = Array{Any}(undef, n_fields);

                         for i=1:n_iter

                              # Set end for the sync loop
                              end_sync_loop = n_workers;
                              if i*n_workers > n_fields
                                   end_sync_loop = n_fields-(i-1)*n_iter;
                              end

                              @sync for j=1:end_sync_loop

                                   ind_par = j+(i-1)*n_workers;

                                   # @async remotecall_fetch
                                   @async par_list[ind_par] =
                                        remotecall_fetch(copy, workers_vec[j], getfield(par_prop, j+(i-1)*n_workers));
                              end
                         end
                    end

                    par = ParSsm(par_list...);

                    # Copy θ_unb and θ_bound in parallel
                    @sync begin
                         @async θ_unb   = remotecall_fetch(copy, workers_vec[1], θ_prop_unb);
                         @async θ_bound = remotecall_fetch(copy, workers_vec[2], θ_prop_bound);
                    end
               end
          end

          # Chains
          chain_θ_unb[:, draw]   = θ_unb;
          chain_θ_bound[:, draw] = θ_bound;


          # ------------------------------------------------------------------------------------------------------------
          # Block 2: Draw alpha (Gibbs step)
          # ------------------------------------------------------------------------------------------------------------

          if draw > burnin && gibbs_step == 1

               # Draw α conditional to the state-space parameters
               α_draw, _                  = kalman_diffuse!(par, 0, 1, 1);
               distr_α[:, :, draw-burnin] = α_draw;
               distr_par[draw-burnin]     = par;

               # Forecast
               density_fcst_draw             = (par.Z*α_draw)';
               density_fcst_draw             = density_fcst_draw[end-h+1:end, :];
               distr_fcst[:, :, draw-burnin] = density_fcst_draw;
          end
     end

     return chain_θ_unb, chain_θ_bound, distr_α, distr_fcst, par, distr_par;
end
