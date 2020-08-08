#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function get_par_bound(θ_unb::Array{Float64,1}, MIN::Array{Float64,1}, MAX::Array{Float64,1}, opt_transf::Array{Int64, 1})

# ----------------------------------------------------------------------------------------------------------------------
# Change domain: from bounded to unbounded domain
# - Transformation: 1 natural logarithm, 2 no transformations, 3 generalized logit
# ----------------------------------------------------------------------------------------------------------------------

     i1 = opt_transf .== 1;
     i2 = opt_transf .== 2;
     i3 = opt_transf .== 3;

     θ_bound     = zeros(size(θ_unb));
     θ_bound[i1] = exp.(θ_unb[i1]) + MIN[i1];
     θ_bound[i2] = θ_unb[i2];
     θ_bound[i3] = (MIN[i3] + MAX[i3] .* exp.(θ_unb[i3])) ./ (1 .+ exp.(θ_unb[i3]));

     return θ_bound;
end
