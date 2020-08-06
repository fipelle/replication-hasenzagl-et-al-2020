#=
Name: tc_main.jl
Description: Execution manager
=#

# ----------------------------------------------------------------------------------------------------------------------
# Execution: run_type == 1
# - Single iteration: it executes the code using the most updated datapoints
# ----------------------------------------------------------------------------------------------------------------------

if run_type == 1

    # ----- Standardize the data -----

    # Mixed-frequency
    if nM > 0 && nQ > 0
         σʸ_monthly   = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM]...);
         σʸ_quarterly = hcat([std(collect(skipmissing(diff(data[3:3:end, i]))), dims=1) for i=nM+1:nM+nQ]...);
         σʸ           = [σʸ_monthly σʸ_quarterly];

    # Mono frequency
    else
         σʸ = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM+nQ]...);
    end

    data = data ./ σʸ; # Standardize the data
    data = [data; missing.*ones(h, nM+nQ)]; # add "missings" for forecast

    # ----- Run the Metropolis-Within-Gibbs -----

    distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, mwg_const, acc_rate, par, par_ind, par_size, distr_par =
         tc_mwg(data, h, nDraws, burnin, mwg_const, σʸ);

    # Save res in jld format
    save(string("./res_", res_name, ".jld"), Dict("distr_α" => distr_α, "distr_fcst" => distr_fcst, "chain_θ_unb" => chain_θ_unb,
                "chain_θ_bound" => chain_θ_bound, "mwg_const" => mwg_const, "acc_rate" => acc_rate, "par" => par,
                "nDraws" => nDraws, "burnin" => burnin, "data" => data, "date" => date, "nM" => nM, "nQ" => nQ,
                "MNEMONIC" => MNEMONIC, "par_ind" => par_ind, "par_size" => par_size, "distr_par" => distr_par, "σʸ" => σʸ));


# ----------------------------------------------------------------------------------------------------------------------
# Execution: run_type == 2
# - Conditional forecasts
# ----------------------------------------------------------------------------------------------------------------------

elseif run_type == 2

    # ----- Load in-sample output -----

    # Load jld output
    res_iis = load(string("./res_", res_name_iis, ".jld"));

    # Minimum output to compute the conditional forecasts
    data      = res_iis["data"];
    date      = res_iis["date"];
    σʸ        = res_iis["σʸ"];
    distr_par = res_iis["distr_par"];
    MNEMONIC  = res_iis["MNEMONIC"];

    # Remove h closing missing values
    data  = data[1:end-(size(data)[1] - size(date)[1]), :];

    # ----- Conditional forecasts -----

    for i=1:length(cond)

        # Size
        observations, n = size(data);

        # Conditioning path
        keys_ith   = collect(keys(cond[i]));
        values_ith = collect(values(cond[i]));
        pos_ith    = vcat([findall(MNEMONIC .== j) for j in keys_ith]...);
        h_ith      = length.(values_ith);

        # Generate conditioning data
        data_ith = copy(data);
        data_ith = [data; missing .* ones(maximum([h; h_ith]), n)];             # add "missings" for conditional forecast
        for j=1:length(pos_ith)
            data_ith[observations+1:observations+h_ith[j], pos_ith[j]] = values_ith[j] ./ σʸ[pos_ith[j]];
        end

        # Pre-allocate memory
        k               = size(distr_par[1].T)[1];
        m, n            = size(data_ith);
        distr_fcst_cond = zeros(m, n, nDraws[2]-burnin[2]);
        distr_α_cond    = zeros(k, m, nDraws[2]-burnin[2]);

        # Run the conditional forecast
        for draw=1:nDraws[2]-burnin[2]
            if draw > 1 && mod(draw, 100) == 0
                print("Conditional forecast $i (out of $(size(cond)[1])) > $draw-th iteration (out of $(nDraws[2]-burnin[2])) \n");
            end

            # Draw
            par_draw                 = distr_par[draw];
            par_draw.y               = permutedims(data_ith);
            α_draw, _                = kalman_diffuse!(par_draw, 0, 1, 1);
            distr_α_cond[:, :, draw] = α_draw;

            # Forecast
            distr_fcst_cond[:, :, draw] = (par_draw.Z * α_draw)' .* σʸ;
        end

        print("\n");

        # Save res in jld format
        save(string("./res_", res_name, "_cond$(i).jld"), Dict("data_ith" => data_ith, "distr_fcst_cond" => distr_fcst_cond,
                                                                "distr_α_cond" => distr_α_cond, "conditional_path" => cond[i]));
    end


# ----------------------------------------------------------------------------------------------------------------------
# Execution: run_type == 3
# -  Out-of-sample: out-of-sample exercise, forecasting period starts after end_presample_vec
# ----------------------------------------------------------------------------------------------------------------------

elseif run_type == 3

    # ----- Initialise -----

    data_full = copy(data);

    # end_presample_vec
    date_vec      = vcat([[Dates.day(date[i]) Dates.month(date[i]) Dates.year(date[i])] for i=1:length(date)]...);
    end_presample = findall(sum(date_vec .== end_presample_vec', dims=2)[:] .== 3)[1];
    end_oos       = size(date_vec)[1];

    if end_presample == end_oos
        error("end_presample_vec is not set correctly");
    end

    oos_forecast = Array{Any}(undef, end_oos-end_presample+1);
    α_array      = Array{Any}(undef, end_oos-end_presample+1);
    σ_array      = Array{Any}(undef, end_oos-end_presample+1);

    # ----- Run the out-of-sample -----

    for t=end_presample:end_oos

        # data and data_cond
        data = data_full[1:t, :];

        # Mixed-frequency
        if nM > 0 && nQ > 0
            σʸ_monthly   = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM]...);
            σʸ_quarterly = hcat([std(collect(skipmissing(diff(data[3:3:end, i]))), dims=1) for i=nM+1:nM+nQ]...);
            σʸ           = [σʸ_monthly σʸ_quarterly];

        # Mono frequency
        else
            σʸ = hcat([std(collect(skipmissing(diff(data[:, i]))), dims=1) for i=1:nM+nQ]...);
        end

        data = data ./ σʸ;
        data = [data; missing .* ones(h, nM+nQ)];

        # Run the Metropolis-Within-Gibbs
        global mwg_const;
        println(mwg_const)
        distr_α, distr_fcst, chain_θ_unb, chain_θ_bound, mwg_const, acc_rate, par, par_ind, par_size, distr_par =
             tc_mwg(data, h, nDraws, burnin, mwg_const, σʸ);

        # Re-attribute standard deviation of the delta to the forecasts
        for draw=1:nDraws[2]-burnin[2]
            distr_fcst[:, :, draw] = distr_fcst[:, :, draw] .* σʸ;
        end

        # Store results
        oos_forecast[t-end_presample+1] = distr_fcst;
        α_array[t-end_presample+1]      = distr_α;
        σ_array[t-end_presample+1]      = σʸ;
    end

    # Save res in jld format
    save(string("./res_", res_name, ".jld"), Dict("distr_α" => distr_α, "distr_fcst" => distr_fcst, "chain_θ_unb" => chain_θ_unb,
                "chain_θ_bound" => chain_θ_bound, "mwg_const" => mwg_const, "acc_rate" => acc_rate, "par" => par,
                "nDraws" => nDraws, "burnin" => burnin, "data" => data, "date" => date, "nM" => nM, "nQ" => nQ,
                "MNEMONIC" => MNEMONIC, "par_ind" => par_ind, "par_size" => par_size, "distr_par" => distr_par,
                "oos_forecast" => oos_forecast, "α_array" => α_array, "σ_array" => σ_array, "data_full" => data_full));
end
