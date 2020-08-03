function get_par_unb(θ_bound::Array{Float64,1}, MIN::Array{Float64,1}, MAX::Array{Float64,1}, opt_transf::Array{Int64, 1})
# ----------------------------------------------------------------------------------------------------------------------
# Change domain: from bounded to unbounded domain
# - Transformation: 1 natural logarithm, 2 no transformations, 3 generalized logit
# ----------------------------------------------------------------------------------------------------------------------

     i1 = opt_transf .== 1;
     i2 = opt_transf .== 2;
     i3 = opt_transf .== 3;

     θ_unb     = zeros(size(θ_bound));
     θ_unb[i1] = log.(θ_bound[i1] - MIN[i1]);
     θ_unb[i2] = θ_bound[i2];
     θ_unb[i3] = log.((θ_bound[i3] - MIN[i3]) ./ (MAX[i3] - θ_bound[i3]));

     return θ_unb;
end
