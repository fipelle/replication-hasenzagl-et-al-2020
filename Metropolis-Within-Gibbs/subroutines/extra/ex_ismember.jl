function ex_ismember(a, b)
     x = sum(a .== b', 2) .== 1;
     return x;
end
