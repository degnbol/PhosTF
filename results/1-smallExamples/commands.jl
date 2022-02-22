#!/usr/bin/env julia
# Run me from one of the folders that has the WT.adj and WP.adj
@src "simulation/simulate";
@src "inference/infer"

nₜ = ncol(loaddlm("WT.adj"; header=true)) - 1
nₚ = ncol(loaddlm("WP.adj"; header=true)) - 1


using .Iterators: product # since @threads does not support nested loops yet
using CSV
using DataFrames
using Logging
incl("src/utilities/AUC")
Logging.disable_logging(Logging.Warn)
Wₜ = loaddlm("WT.adj"; header=true)
# make sure to enforce that it is Int which indicates weight presence and sign as opposed to Float that indicates weight magnitude.
Wₚ = loaddlm("WP.adj", Int; header=true)
Wₜ = Matrix(Wₜ[:, 2:end])
Wₚ = Matrix(Wₚ[:, 2:end])

open("aucs.tsv", "w") do io
    header = ["cancel", "mean_k", "vec", "norm", "lam", "AUC_T", "AUC_P"]
    write(io, join(header, '\t') * '\n')
    for (i, _cancel, _mean_k, _vec, _norm, _lam) in collect(product(1:5, [false, true], [true, false], [true, false], 0:2, 0:4))
        hyperD = Dict{String,Real}("cancel"=>_cancel, "mean_k"=>_mean_k, "vec"=>_vec, "norm"=>_norm, "lam"=>_lam)
        for _λB in [0., 0.01, 0.1, 0.5, 1.]
            for _λW in [0., 0.01, 0.1, 0.5, 1.]
                net = GeneRegulation.Network(Wₜ, Wₚ; hyper=hyperD)
                logFC = ODEs.@domainerror(ODEs.logFC(net))
                if measurements !== nothing
                    savedlm("sim_logFC.tsv", measurements; colnames=net.names[1:net.nₜ+net.nₚ], rownames=net.names)
                    infer("sim_logFC.tsv", r"T.*", r"P.*"; epochs=1000, lambda_Bstar=_λB, lambda_absW=_λW)
                    WT_infer = loaddlm("WT_infer.tsv"; header=true)[!, 2:end]
                    WP_infer = loaddlm("WP_infer.tsv"; header=true)[!, 2:end]
                    hyperD["AUC_T"] = AUC(Matrix(WT_infer), Matrix(WT))
                    hyperD["AUC_P"] = AUC(Matrix(WP_infer), Matrix(WP))
                else
                    hyperD["AUC_T"] = NaN
                    hyperD["AUC_P"] = NaN
                end
                write(io, join([hyperD[k] for k in header], '\t') * '\n')
            end
        end
    end
end

# 
# plot_timeseries("sim.tsv", nₜ, nₚ; o="sim.pdf")
# for mut_id in 1:nₜ+nₚ
#     timeseries("net.bson", mut_id; u0="sim.tsv")
#     plot_timeseries("sim_mut$mut_id.tsv", nₜ, nₚ; o="sim_mut$mut_id.pdf")
# end
# steadystate("net.bson")
# for mut_id in 1:nₜ+nₚ
#     steadystate("net.bson", mut_id; u0="sim.tsv")
# end
# logFC("net.bson", "sim_logFC.tsv")
# 
# @src "inference/infer"
# infer("sim_logFC.tsv", r"T.*", r"P.*"; epochs=10000, lambda_Bstar=.1, lambda_absW=1.)
# 
# @src "visualization/W2graph"
# xgmml("net.bson"; o="net.xgmml")
# xgmml("net.bson"; o="steady.xgmml", X="steady.tsv")
# for mut_id in 1:nₜ+nₚ
#     xgmml("net.bson"; o="steady_mut$mut_id.xgmml", X="steady_mut$mut_id.tsv", highlight=mut_id)
# end
# xgmml("net.bson"; o="ko.xgmml", X="sim_logFC.tsv")
# xgmml("WT_infer.tsv", "WP_infer.tsv"; o="infer.xgmml")
# 
# 

