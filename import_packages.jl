#=
This file is part of the replication code for: Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation.
Please cite the paper if you are using any part of the code for academic work (including, but not limited to, conference and peer-reviewed papers).
=#

import Pkg;

# Packages used in the estimation
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("FileIO")
Pkg.add("JLD")
Pkg.add("XLSX")

# Packages used to make the charts
Pkg.add("Colors")
Pkg.add("CSV")
Pkg.add("ORCA")
Pkg.add("KernelDensity")
Pkg.add("PlotlyJS")
