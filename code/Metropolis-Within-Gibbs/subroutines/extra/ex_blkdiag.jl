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
