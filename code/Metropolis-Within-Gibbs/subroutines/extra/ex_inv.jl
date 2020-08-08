#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function ex_inv(x)

# ----------------------------------------------------------------------------------------------------------------------
# Custom inv function
# ----------------------------------------------------------------------------------------------------------------------

    xInv = x\eye(size(x)[1], size(x)[2]);

    return xInv;
end
