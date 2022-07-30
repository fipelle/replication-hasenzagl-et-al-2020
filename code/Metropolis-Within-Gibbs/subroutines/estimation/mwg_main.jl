#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function mwg_main(par::ParSsm, h::Int64, nDraws::Array{Int64, 1}, burnin::Array{Int64, 1},
                  mwg_const::Array{Float64, 1}, par_ind::BoolParSsm, t=0::Int64, end_oos=0::Int64)

# ----------------------------------------------------------------------------------------------------------------------
# Metropolis-Within-Gibbs algorithm: mainframe
# ----------------------------------------------------------------------------------------------------------------------

     # -----------------------------------------------------------------------------------------------------------------
     # Error management
     # -----------------------------------------------------------------------------------------------------------------

     for j in axes(par_ind.Z, 2), i in axes(par_ind.Z, 1)
          if (par_ind.Z[i,j] + par_ind.Z_plus[i,j] + par_ind.Z_minus[i,j] + par_ind.Z_bounded[i,j]) > 1
               error("Each measurement equation coefficient can be either unrestricted or constrained - not both!");
          end
     end

     
     # -----------------------------------------------------------------------------------------------------------------
     # Initialise
     # -----------------------------------------------------------------------------------------------------------------

     # Dimensions
     k    = size(par.T)[1];
     n, m = size(par.y);

     # Target acceptance rate
     acc_target = 30;

     # No. parameters to be estimated
     par_size = SizeParSsm(sum(par_ind.d),
                                 sum(sum(par_ind.Z)),
                                 sum(sum(par_ind.Z_plus)),
                                 sum(sum(par_ind.Z_minus)),
                                 sum(sum(par_ind.Z_bounded)),
                                 sum(sum(par_ind.R)),
                                 sum(par_ind.c),
                                 sum(sum(par_ind.T)),
                                 sum(sum(par_ind.Q)),
                                 sum(sum(par_ind.λ)),
                                 sum(sum(par_ind.ρ)),
                                 sum(par_ind.d) + sum(sum(par_ind.Z)) + sum(sum(par_ind.Z_plus)) +
                                    sum(sum(par_ind.Z_minus)) + sum(sum(par_ind.Z_bounded)) + sum(sum(par_ind.R)) +
                                    sum(par_ind.c) + sum(sum(par_ind.T)) + sum(sum(par_ind.Q)) +
                                    sum(sum(par_ind.λ)) + sum(sum(par_ind.ρ)));


     # -----------------------------------------------------------------------------------------------------------------
     # Lower and upper bounds
     # -----------------------------------------------------------------------------------------------------------------

     # Numerical constants
     xi                = 1e-3; # it must be small
     MIN_var           = 0;
     MIN_coeff         = -Inf;
     MIN_coeff_plus    = 0;
     MIN_coeff_minus   = -Inf;
     MIN_coeff_bounded = -1;
     MIN_λ             = xi;
     MIN_ρ             = xi;
     MAX_var           = Inf;
     MAX_coeff         = Inf;
     MAX_coeff_plus    = Inf;
     MAX_coeff_minus   = 0;
     MAX_coeff_bounded = 1;
     MAX_λ             = pi;
     MAX_ρ             = 0.97;

     # MIN
     MIN = [MIN_var*ones(par_size.R);
            MIN_coeff*ones(par_size.d + par_size.Z);
            MIN_coeff_plus*ones(par_size.Z_plus);
            MIN_coeff_minus*ones(par_size.Z_minus);
            MIN_coeff_bounded*ones(par_size.Z_bounded);
            MIN_var*ones(par_size.Q);
            MIN_coeff*ones(par_size.c + par_size.T);
            MIN_λ*ones(par_size.λ);
            MIN_ρ*ones(par_size.ρ)];

     # MAX
     MAX = [MAX_var*ones(par_size.R);
            MAX_coeff*ones(par_size.d + par_size.Z);
            MAX_coeff_plus*ones(par_size.Z_plus);
            MAX_coeff_minus*ones(par_size.Z_minus);
            MAX_coeff_bounded*ones(par_size.Z_bounded);
            MAX_var*ones(par_size.Q);
            MAX_coeff*ones(par_size.c + par_size.T);
            MAX_λ*ones(par_size.λ);
            MAX_ρ*ones(par_size.ρ)];


     # -----------------------------------------------------------------------------------------------------------------
     # Set prior parameters
     # -----------------------------------------------------------------------------------------------------------------

     # Define prior distribution objects
     prior_opt = PriorOpt(Normal(0, 1/xi),
                          TruncatedNormal(0, 1/xi, MIN_coeff_plus, MAX_coeff_plus),
                          TruncatedNormal(0, 1/xi, MIN_coeff_minus, MAX_coeff_minus),
                          TruncatedNormal(0, 1/xi, MIN_coeff_bounded, MAX_coeff_bounded),
                          InverseGamma(3, 1),
                          par_size.λ*logpdf.(Uniform(MIN_λ, MAX_λ), MIN_λ),
                          par_size.ρ*logpdf.(Uniform(MIN_ρ, MAX_ρ), MIN_ρ));

     # Transformations: 1 natural logarithm, 2 no transformations, 3 generalized logit
     opt_transf = convert(Array{Int64, 1}, [1*ones(par_size.R);
                                            2*ones(par_size.d + par_size.Z);
                                            1*ones(par_size.Z_plus);
                                            0*ones(par_size.Z_minus);
                                            3*ones(par_size.Z_bounded);
                                            1*ones(par_size.Q);
                                            2*ones(par_size.c + par_size.T);
                                            3*ones(par_size.λ);
                                            3*ones(par_size.ρ)]);


     # -----------------------------------------------------------------------------------------------------------------
     # Initialisation of the Metropolis-Within-Gibbs algorithm
     # -----------------------------------------------------------------------------------------------------------------

     # Initial guess for θ_ini_bound
     θ_ini_bound = [ones(par_size.R);
                    zeros(par_size.d);
                    ones(par_size.Z + par_size.Z_plus);
                    -ones(par_size.Z_minus);
                    (MAX_coeff_bounded-MIN_coeff_bounded)*ones(par_size.Z_bounded);
                    ones(par_size.Q);
                    zeros(par_size.c);
                    ones(par_size.T);
                    (2*pi/32)*ones(par_size.λ);
                    0.5*ones(par_size.ρ)];

     # Run init
     θ_ini_unb = get_par_unb(θ_ini_bound, MIN, MAX, opt_transf);

     algorithm_name = "Initialisation";
     acc_rate       = 0.0;

     chain_θ_unb_init = Array{Any}(undef, 1);
     Σ_init           = Array{Float64}(Matrix(I, length(θ_ini_unb), length(θ_ini_unb)));

     while acc_rate < 25.0 || acc_rate > 35.0

          chain_θ_unb_init, _, _, _, _, _ =
               mwg_run(θ_ini_unb, par, h, par_ind, par_size, prior_opt, MIN, MAX, opt_transf, nDraws[1], burnin[1],
                    mwg_const[1], algorithm_name, t, end_oos, Σ_init, 0);

          acc_rate = get_progress(chain_θ_unb_init, burnin[1]+1, nDraws[1], 0);

          if acc_rate < 25.0 || acc_rate > 35.0
               delta_acc    = acc_rate-acc_target;
               mwg_const[1] = mwg_const[1]*exp(delta_acc/100);
          end
     end

     θ_unb_init_clean = chain_θ_unb_init[:, burnin[1]+1:end]';
     θ_start_unb      = median(θ_unb_init_clean, dims=1)[:];
     Σ_mwg            = cov(θ_unb_init_clean, dims=1);

     # -----------------------------------------------------------------------------------------------------------------
     # Metropolis-Within-Gibbs algorithm
     # -----------------------------------------------------------------------------------------------------------------

     algorithm_name = "Metropolis-Within-Gibbs";
     acc_rate       = 0.0;

     chain_θ_unb   = Array{Any}(undef, 1);
     chain_θ_bound = Array{Any}(undef, 1);
     distr_α       = Array{Any}(undef, 1);
     distr_fcst    = Array{Any}(undef, 1);
     distr_par     = Array{Any}(undef, 1);

     while acc_rate < 25.0 || acc_rate > 35.0

          chain_θ_unb, chain_θ_bound, distr_α, distr_fcst, par, distr_par =
               mwg_run(θ_start_unb, par, h, par_ind, par_size, prior_opt, MIN, MAX, opt_transf, nDraws[2], burnin[2],
                    mwg_const[2], algorithm_name, t, end_oos, Σ_mwg, 1);

          acc_rate = get_progress(chain_θ_unb, burnin[2]+1, nDraws[2], 0);

          if acc_rate < 25.0 || acc_rate > 35.0
               delta_acc    = acc_rate-acc_target;
               mwg_const[2] = mwg_const[2]*exp(delta_acc/100);
          end
     end

     return distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, mwg_const, acc_rate, par, par_size, distr_par;
end
