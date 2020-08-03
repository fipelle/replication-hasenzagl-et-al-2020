## A Model of the Fed's View on Inflation

This repository contains the source code for replicating the results in the paper:

Hasenzagl, T., Pellegrino, F., Reichlin, L., & Ricco, G. (2020). A Model of the Fed's View on Inflation

If you have any questions, comments, or suggestions please create a new issue.

## Code structure
The main directory is organized as follows:

* *code_main*: Contains all of the julia code necessary for replication.
    + The *Metropolis-Within-Gibbs* subdirectory contains the code for the Metropolis-Within-Gibbs algorithm.
* *data*: Contains the data used in the estimation. The data is saved in ".csv" and ".xlsx" files. 
* *csv_output*: Contains the generated csv output files.   
* *img*: Contains the generated figures.

## Running the code

The code was written and run in Julia 1.0: https://julialang.org/

The code uses a number of Julia packages. All necessary packages can be installed using the `import_packages.jl` script. To do so, start julia and use the following command at the julia REPL prompt: `julia> include("import_packages.jl")`

The main file is `user_main.jl`. This script runs the

* in-sample estimation, by setting `run_type=1` in `user_main.jl`.
* conditional forecasting exercise, by setting `run_type=2` in `user_main.jl` and specifying the start date of the forecasting exercise, and the conditioning variables and time periods. 
* out-of-sample forecasting exercise, by setting `run_type=3` and specifying the start date of the forecasting exercise.

To run the script start julia and use the following command at the julia REPL prompt: ``julia> include("user_main.jl")``

## Figures and Tables

The figures and tables are created in two Jupyter (https://jupyter.org/) notebooks:

* `iis_charts.ipynb`: creates all figures relating to the in-sample estimation
* `oos_charts.ipynb`: creates all figures relating to the out-of-sample forecasting exercise and the RMSE of the trend-cycle model relative to the RMSE of a random walk with drift. 
