function get_progress(chain::Array{Float64, 2}, istart::Int64, iend::Int64, print_progress=1::Int64, acc_str=""::String)
# ----------------------------------------------------------------------------------------------------------------------
# Get acceptance rate and print current status of the estimation (if enabled)
# ----------------------------------------------------------------------------------------------------------------------

     chain_iter = chain[:, istart:iend];
     acc_rate         = 100*mean(chain_iter[1, 2:end] .!= chain_iter[1, 1:end-1]);
     acc_rate         = round(acc_rate, digits=2);

     if print_progress == 1
          if acc_str != ""
               print("- $acc_str: $acc_rate%\n");
          else
               print("- Acceptance rate: $acc_rate%\n");
          end
     end

     return acc_rate;
end
