function get_mwg_jump(Σ::Array{Float64, 2}, mwg_const::Float64, nDraws_worker::Int64)
# ----------------------------------------------------------------------------------------------------------------------
# Get the Metropolis-Within-Gibbs jumps for the proposal distribution
# ----------------------------------------------------------------------------------------------------------------------

     len_Σ              = size(Σ)[1];
     jumps_obj          = MvNormal(zeros(len_Σ), (mwg_const^2)*Σ);
     θ_unb_jumps_worker = zeros(len_Σ, nDraws_worker);

     for draw=1:nDraws_worker
          θ_unb_jumps_worker[:, draw] = rand(jumps_obj);
     end

     return θ_unb_jumps_worker;
end
