#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function kalman_diffuse!(par::ParSsm, do_loglik=0::Int64, do_smoother=0::Int64, do_sim_smoother=0::Int64, F_tol=1e-8::Float64)

# ----------------------------------------------------------------------------------------------------------------------
# Diffuse Kalman Filter, Smoother and Simulation Smoother (univariate approach for multivariate models).
# The function update the loglikelihood in par if do_loglik is 1.
#
# The model is specified as follows:
#
# y_{t} = d + Z*a_{t}   + e_{t}, where e_{t} ~ N(0, R)
# α_{t} = c + T*α_{t-1} + u_{t}, where u_{t} ~ N(0, Q)
#
# ----------------------------------------------------------------------------------------------------------------------
#
# Kalman routines for non-stationary systems with time-varying system matrices and missing data. The Durbin and
# Koopman's simulation smoother includes the Jarocinski (2015) modifications. The univariate approach for multivariate
# models introduced by Anderson and Moore (1979) and extended by Durbin and Koopman (2012) is used to improve the
# estimation speed.
#
# See: Harvey Koopman and Shephard (2004) or Durbin and Koopman (2012)
#
# ----------------------------------------------------------------------------------------------------------------------
# Input arguments:
#
# y       observed data                                              (n x m)
# d       constant of the observation equation                       (n x 1)
# Z       coefficients of the observation equation                   (n x k)
# R       variance-covariance matrix of the observation equation     (n x n)
# c       constant of the transition equation                        (k x 1)
# T       coefficients of the transition equation                    (k x k)
# Q       variance-covariance matrix of the transition equation      (k x k)
# α¹      initial value of the states                                (k x 1)
# P¹      initial variance-covariance matrix of the states           (k x k)
# P̄¹      diagonal matrix: 1 for non-stationary states (0 otherwise) (k x k)
# F_tol   small value used to determine whether F is zero or not     (1 x 1)
#
# ----------------------------------------------------------------------------------------------------------------------
#
# Optional input arguments:
#
# loglik_bool    enable the estimation of the loglikelihood        (1 x 1)
# sim_smoother   enable the simulation smoother                    (1 x 1)
#
# ----------------------------------------------------------------------------------------------------------------------
#
# Output
# α̃              smoothed state vector                             (k x m)
# Pˢ             smoothed variance-covariance matrix of the states (k x k)
# loglik         loglikelihood (if enabled)                        (1 x 1)
# ----------------------------------------------------------------------------------------------------------------------

     # -----------------------------------------------------------------------------------------------------------------
     # Initialisation
     # -----------------------------------------------------------------------------------------------------------------

     # Dimensions
     k    = size(par.T)[1]
     n, m = size(par.y);

     # Kalman filter output
     ν  = zeros(m, n);
     K  = zeros(k, n, m);
     K̄  = zeros(k, n, m);
     F  = zeros(m, n);
     F̄  = zeros(m, n);
     α̂ᶠ = zeros(k, m);
     Pᶠ = zeros(k, k, m);
     P̄ᶠ = zeros(k, k, m);

     # Kalman smoother output
     α̂ˢ = zeros(k, m);
     Pˢ = zeros(k, k, m);
     P̄ˢ = zeros(k, k, m);

     # Loglikelihood
     if do_loglik == 1
          par.loglik = 0;
     end

     if do_sim_smoother == 1
          do_smoother = 1;
     end

     # Initial variance
     par.P¹[par.P̄¹ .== 1] .= 0;


     # -----------------------------------------------------------------------------------------------------------------
     # Simulation smoother (if enabled)
     # -----------------------------------------------------------------------------------------------------------------

     if do_sim_smoother == 1

          # Draw: y⁺ and α⁺

          y⁺ = zeros(n, m);
          α⁺ = zeros(k, m);

          # ------------------------------------------------------------------------------------------------------------
          # Durbin, Koopman 2002, Algorithm 2.
          #
          # Generate y⁺ and α⁺ - y and α are drawn from their unconditional distributions using the demeaned model,
          # i.e. zero initial state and zero constant terms.
          #
          # Note: Jarocinski (2015) modification
          # ------------------------------------------------------------------------------------------------------------

          for t=1:m
               if t == 1
                    # draw draw the first state with α1=0
                    α⁺[:, 1] = get_random_disturbance(k, par.P¹);
               else
                    α⁺[:, t] = par.T*α⁺[:, t-1] + get_random_disturbance(k, par.Q);
               end

               y⁺[:, t] = par.Z*α⁺[:, t] + get_random_disturbance(n, par.R);
          end

          # Define antithetic variable yᵃ
          yᵃ = par.y - y⁺;

     else
          yᵃ = par.y;
     end


     # -----------------------------------------------------------------------------------------------------------------
     # In order to use the univariate approach we must deal with correlated measurement errors (if any)
     # -----------------------------------------------------------------------------------------------------------------

     # 1. Check if there are correlated measurement errors
     R_ind               = par.R .!= 0;
     R_ind[Matrix(I,n,n) .== 1] .= false;

     # 2. If there are correlated measurement errors, apply the following Shur decomposition R = S⋅H⋅S', where H is
     #    diagonal and holds the eigenvalues of R, while S is orthogonal. Pre-multiply the observation equation by S',
     #    so that the variance of the errors will be diagonal and equal to H
     if sum(sum(R_ind)) > 0
          H, S, _ = schur(copy(par.R));
          yᵃ = S'*yᵃ;
          d̃  = S'*copy(par.d);
          Z̃  = S'*copy(par.Z);
          R̃  = H;
     else
          d̃ = par.d;
          Z̃ = par.Z;
          R̃ = par.R;
     end


     # -----------------------------------------------------------------------------------------------------------------
     # Kalman Filter - Univariate (or sequential) filtering
     # -----------------------------------------------------------------------------------------------------------------

     α̂ᶠᵢ = par.α¹;
     Pᶠᵢ = par.P¹;
     P̄ᶠᵢ = par.P̄¹;

     α̂ᶠ[:, 1]    = α̂ᶠᵢ;
     Pᶠ[:, :, 1] = Pᶠᵢ;
     P̄ᶠ[:, :, 1] = P̄ᶠᵢ;

     # Loop over t
     for t=1:m

          # Get last estimates for α and P
          α̂ᶠᵢ = α̂ᶠ[:, t];
          Pᶠᵢ = Pᶠ[:, :, t];
          P̄ᶠᵢ = P̄ᶠ[:, :, t];

          # Handle missing data
          no_na     = ismissing.(yᵃ[:, t]) .== false;
          ind_no_na = findall(no_na);
          yᵃᵗ       = convert(Array{Float64, 1}, yᵃ[no_na, t]);
          d̃ᵗ        = d̃[no_na, :];
          Z̃ᵗ        = Z̃[no_na, :];
          R̃ᵗ        = R̃[no_na, no_na];

          # Loop over not missing data
          for i=1:sum(no_na)

               # Kalman filter i-th output
               Fᵗᵢ = (Z̃ᵗ[i, :]'*Pᶠᵢ*Z̃ᵗ[i, :] + R̃ᵗ[i, i])[1];
               F̄ᵗᵢ = (Z̃ᵗ[i, :]'*P̄ᶠᵢ*Z̃ᵗ[i, :])[1];

               if F̄ᵗᵢ > F_tol || Fᵗᵢ > F_tol
                    Kᵗᵢ = Pᶠᵢ*Z̃ᵗ[i, :];
                    K̄ᵗᵢ = P̄ᶠᵢ*Z̃ᵗ[i, :];
                    νᵗᵢ = (yᵃᵗ[i, :] .- d̃ᵗ[i, :] .- Z̃ᵗ[i, :]'*α̂ᶠᵢ)[1];


                    # --------------------------------------------------------------------------------------------------
                    # Case 1: F̄ᵗᵢ greater than zero
                    # --------------------------------------------------------------------------------------------------

                    if F̄ᵗᵢ > F_tol

                         if t < m
                              # Update α̂ᶠᵢ, Pᶠᵢ and P̄ᶠᵢ
                              α̂ᶠᵢ = α̂ᶠᵢ + K̄ᵗᵢ/F̄ᵗᵢ*νᵗᵢ;
                              Pᶠᵢ = Pᶠᵢ + K̄ᵗᵢ*K̄ᵗᵢ'*Fᵗᵢ/(F̄ᵗᵢ^2) - (Kᵗᵢ*K̄ᵗᵢ'+K̄ᵗᵢ*Kᵗᵢ')/F̄ᵗᵢ;
                              P̄ᶠᵢ = P̄ᶠᵢ - K̄ᵗᵢ/F̄ᵗᵢ*K̄ᵗᵢ';
                         end

                         # Estimate loglikelihood  (if enabled)
                         if do_loglik == 1
                              par.loglik = par.loglik -0.5*(log(2*pi) + log(F̄ᵗᵢ));
                         end


                    # --------------------------------------------------------------------------------------------------
                    # Case 2: F̄ᵗᵢ is equal zero
                    # --------------------------------------------------------------------------------------------------

                    else

                         if t < m
                              # Update α̂ᶠᵢ, Pᶠᵢ and P̄ᶠᵢ (the latter doesn't change in this case)
                              α̂ᶠᵢ = α̂ᶠᵢ + Kᵗᵢ/Fᵗᵢ*νᵗᵢ;
                              Pᶠᵢ = Pᶠᵢ - Kᵗᵢ/Fᵗᵢ*Kᵗᵢ';
                         end

                         # Estimate loglikelihood  (if enabled)
                         if do_loglik == 1
                              par.loglik = par.loglik -0.5*(log(2*pi) + log(Fᵗᵢ) + (νᵗᵢ^2)/Fᵗᵢ);
                         end
                    end


                    # Store Kalman filter output (1/2): i-th values
                    if do_smoother == 1
                         ν[t, ind_no_na[i]]    = νᵗᵢ;
                         K[:, ind_no_na[i], t] = Kᵗᵢ;
                         K̄[:, ind_no_na[i], t] = K̄ᵗᵢ;
                         F[t, ind_no_na[i]]    = Fᵗᵢ;
                         F̄[t, ind_no_na[i]]    = F̄ᵗᵢ;
                    end
               end

               if t < m
                    # Store Kalman filter output (2/2): t-th values
                    α̂ᶠ[:, t+1]    = par.c + par.T*α̂ᶠᵢ;
                    Pᶠ[:, :, t+1] = par.T*Pᶠᵢ*par.T' + par.Q;
                    P̄ᶠ[:, :, t+1] = par.T*P̄ᶠᵢ*par.T';
               end
          end

          if sum(no_na) == 0 && t < m

               # Store Kalman filter predictions
               α̂ᶠ[:, t+1]    = par.c + par.T*α̂ᶠᵢ;
               Pᶠ[:, :, t+1] = par.T*Pᶠᵢ*par.T' + par.Q;
               P̄ᶠ[:, :, t+1] = par.T*P̄ᶠᵢ*par.T';
          end
     end


     # -----------------------------------------------------------------------------------------------------------------
     # Kalman Smoother - Univariate (or sequential) smoother
     # -----------------------------------------------------------------------------------------------------------------

     if do_smoother == 1

          rᵢ = zeros(2*k);
          Nᵢ = zeros(2*k, 2*k);
          TT = kron(Matrix(I,2,2), par.T);

          # Loop over t
          for t=m:-1:1

               if t < m
                    rᵢ = TT'*rᵢ;
                    Nᵢ = TT'*Nᵢ*TT;
               end

               # Handle missing data
               no_na     = ismissing.(yᵃ[:, t]) .== false;
               ind_no_na = findall(no_na);
               yᵃᵗ       = convert(Array{Float64, 1}, yᵃ[no_na, t]);
               Z̃ᵗ        = Z̃[no_na, :];

               # Loop over not missing data
               for i=sum(no_na):-1:1

                    # Retrieve Kalman Filter output (1/2)
                    Fᵗᵢ = F[t, ind_no_na[i]];
                    F̄ᵗᵢ = F̄[t, ind_no_na[i]];

                    if F̄ᵗᵢ > F_tol || Fᵗᵢ > F_tol

                         # Retrieve Kalman Filter output (2/2)
                         νᵗᵢ = ν[t, ind_no_na[i]];
                         Kᵗᵢ = K[:, ind_no_na[i], t];
                         K̄ᵗᵢ = K̄[:, ind_no_na[i], t];


                         # ---------------------------------------------------------------------------------------------
                         # Case 1: F̄ᵗᵢ greater than zero
                         # ---------------------------------------------------------------------------------------------

                         if F̄ᵗᵢ > F_tol

                              # Estimate Lᵗᵢ, L̄ᵗᵢ and Mᵗᵢ using the Kalman Filter output
                              Lᵗᵢ = (K̄ᵗᵢ*Fᵗᵢ/F̄ᵗᵢ-Kᵗᵢ)*Z̃ᵗ[i, :]'/F̄ᵗᵢ;
                              L̄ᵗᵢ = Matrix(I,k,k) - K̄ᵗᵢ*Z̃ᵗ[i, :]'/F̄ᵗᵢ;
                              Mᵗᵢ = [L̄ᵗᵢ Lᵗᵢ; zeros(size(L̄ᵗᵢ)) L̄ᵗᵢ];

                              # Estimate rᵢ and Nᵢ
                              temp¹ = Z̃ᵗ[i, :]/F̄ᵗᵢ;
                              temp² = temp¹*Z̃ᵗ[i, :]';
                              rᵢ    = [zeros(k); temp¹*νᵗᵢ] + Mᵗᵢ'*rᵢ;
                              Nᵢ    = [zeros(k, k) temp²; temp² temp²*Fᵗᵢ/F̄ᵗᵢ] + Mᵗᵢ'*Nᵢ*Mᵗᵢ;


                         # ---------------------------------------------------------------------------------------------
                         # Case 2: F̄ᵗᵢ is equal zero
                         # ---------------------------------------------------------------------------------------------

                         else
                              # Estimate Lᵗᵢ, L̄ᵗᵢ and LMᵗᵢ using the Kalman Filter output
                              Lᵗᵢ = Matrix(I,k,k) - Kᵗᵢ*Z̃ᵗ[i, :]'/Fᵗᵢ;
                              Mᵗᵢ = [Lᵗᵢ zeros(size(Lᵗᵢ)); zeros(size(Lᵗᵢ)) Lᵗᵢ];

                              # Estimate rᵢ and Nᵢ
                              temp¹ = Z̃ᵗ[i, :]/Fᵗᵢ;
                              rᵢ    = [temp¹*νᵗᵢ; zeros(k)] + Mᵗᵢ'*rᵢ;
                              Nᵢ    = [temp¹*Z̃ᵗ[i, :]' zeros(k, k); zeros(k, 2*k)] + Mᵗᵢ'*Nᵢ*Mᵗᵢ;
                         end
                    end
               end

               # Smoothed states and covariance
               P̃ᶠ          = [Pᶠ[:, :, t] P̄ᶠ[:, :, t]];
               α̂ˢ[:, t]    = α̂ᶠ[:, t] + P̃ᶠ*rᵢ;
               Pˢ[:, :, t] = Pᶠ[:, :, t] - P̃ᶠ*Nᵢ*P̃ᶠ';
          end
     end


     # -----------------------------------------------------------------------------------------------------------------
     # Return output
     # -----------------------------------------------------------------------------------------------------------------

     if do_sim_smoother == 1
          α = α̂ˢ + α⁺;
          P = Pˢ;

     elseif do_smoother == 1
          α = α̂ˢ;
          P = Pˢ;

     else
          α = α̂ᶠ;
          P = Pᶠ;
     end

     return α, P;
end
