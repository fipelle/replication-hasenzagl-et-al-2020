#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function tc_mwg(y, h, nDraws, burnin, mwg_const, σʸ)

# ----------------------------------------------------------------------------------------------------------------------
# Define the basic structure of the state-space parameters
# ----------------------------------------------------------------------------------------------------------------------

     n = size(y)[2];


     # -----------------------------------------------------------------------------------------------------------------
     # Observation equations
     # -----------------------------------------------------------------------------------------------------------------

     d   = zeros(n);
     Z   = [ones(n) [0; ones(n-1)] [zeros(n-2); ones(2)] [zeros(3); ones(n-3)] [zeros(4); ones(n-4)] [zeros(4); 1 ./ σʸ[end-3:end]]];
     Z1a = kron(Matrix(I, 4, 4), [1, 0, 1])';                                                                                     # idio C, idio C+, idio trend
     Z1b = kron(Matrix(I, 2, 2), [1, 0, 1])';                                                                                     # idio C, idio C+, idio trend
     Z2  = kron(Matrix(I, 2, 2), [1, 0])';                                                                                                    # idio C, idio C+

     Z  = [Z ex_blkdiag(Z1a, Z2, Z1b)];
     R  = zeros(n, n);                                                                                                   # irregular components

     # Indeces for observation equations
     d_ind = d .!= 0;
     Z_ind = zeros(size(Z)) .!= 0;
     R_ind = R .!= 0;

     # Projections
     Z_ind[2:end, [1,2]] .= true; # All        -> PC cycle, t and t-1
     Z_ind[7:end, 3]     .= true; # Expect.    -> PC cycle, t-2
     Z_ind[5:end, [4,5]] .= true; # Prices     -> EP cycle, t and t-1


     # -----------------------------------------------------------------------------------------------------------------
     # Transition equations
     # -----------------------------------------------------------------------------------------------------------------

     c              = zeros(size(Z)[2]);
     ind_trends     = [9; 12]; # GDP, EMPL
     c[ind_trends] .= 1;  # random walk drift

     T_c     = convert(Array{Float64, 2}, [1 0; 0 0]);
     T_ct    = convert(Array{Float64, 2}, [1 0 0; 0 0 0; 0 0 1]);
     T_c_ext = convert(Array{Float64, 2}, [1 0 0; 0 0 0; 0 1 0]);
     Q_c_ext = convert(Array{Float64, 2}, [1 0 0; 0 0 0; 0 0 0]);

     T = cat(dims=[1,2], T_c_ext, [T_ct for i=1:5]..., [T_c for i=1:2]..., [T_ct for i=1:2]...);
     Q = cat(dims=[1,2], Q_c_ext, [T_ct for i=1:5]..., [T_c for i=1:2]..., [T_ct for i=1:2]...);

     # Indeces for transition equations
     c_ind = c .!= 0;
     T_ind = zeros(size(T)) .== 1;
     Q_ind = Q .== 1;

     # Initial conditions for the non-stationary states
     P̄_c   = convert(Array{Float64, 2}, [0 0; 0 0]);
     P̄_ct  = convert(Array{Float64, 2}, [0 0 0; 0 0 0; 0 0 1]);
     P̄¹    = cat(dims=[1,2], zeros(3,3), [P̄_ct for i=1:5]..., [P̄_c for i=1:2]..., [P̄_ct for i=1:2]...);

     # Initial conditions
     α¹ = zeros(size(c));
     P¹       = zeros(size(P̄¹));

     # Trigonometric states
     λ_c   = convert(Array{Float64, 1}, [1; 0]);
     λ_ct  = convert(Array{Float64, 1}, [1; 0; 0]);
     λ     = vcat([1; 0; 0], [λ_ct for i=1:5]..., [λ_c for i=1:2]..., [λ_ct for i=1:2]...);
     ρ     = copy(λ);
     λ_ind = λ .!= 0;
     ρ_ind = copy(λ_ind);


     # -----------------------------------------------------------------------------------------------------------------
     # Metropolis-Within-Gibbs
     # -----------------------------------------------------------------------------------------------------------------

     par_ind = BoolParSsm(d_ind, Z_ind, R_ind, c_ind, T_ind, Q_ind, λ_ind, ρ_ind);
     par     = ParSsm(permutedims(y), d, Z, R, c, T, Q, α¹, P¹, P̄¹, λ, ρ, 0.0, 0.0, 0.0);

     distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, mwg_const, acc_rate, par, par_size, distr_par =
          mwg_main(par, h, nDraws, burnin, mwg_const, par_ind);

     return distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, mwg_const, acc_rate, par, par_ind, par_size, distr_par;
end
