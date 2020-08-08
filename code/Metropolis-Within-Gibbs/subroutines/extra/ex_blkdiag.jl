#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

function ex_blkdiag(x, y, z=[])

# ----------------------------------------------------------------------------------------------------------------------
# Custom blkdiag function
# ----------------------------------------------------------------------------------------------------------------------

    # Initialise
    rx, cx = size(x);
    ry, cy = size(y);

    if size(z)[1] == 0
        rz = 0;
        cz = 0;
    else
        rz, cz = size(z);
    end

    # Blkdiag
    out                                     = zeros(rx+ry+rz, cx+cy+cz);
    out[1:rx, 1:cx]                         = x;
    out[rx+1:rx+ry, cx+1:cx+cy]             = y;
    out[rx+ry+1:rx+ry+rz, cx+cy+1:cx+cy+cz] = z;

    return out;
end
