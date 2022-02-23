#!/usr/bin/env julia
# Run me from one of the folders that has the WT.adj and WP.adj
using IterTools: product
using Printf
@src "simulation/simulate";

nₜ = ncol(loaddlm("WT.adj"; header=true)) - 1
nₚ = ncol(loaddlm("WP.adj"; header=true)) - 1

ideal = loaddlm("ideal_logFC.tsv"; header=true)[:, 2:end] |> Matrix |> vec

thres = 0.1


open("idealScore.tsv", "w") do io
    header = ["cancel", "mean_k", "vec", "norm", "lam", "score"]
    write(io, join(header, '\t') * '\n')
    toIter = collect(product(1:5, [false, true], [true, false], [true, false], 0:2, 0:4))
    progress = 0
    for (i, _cancel, _mean_k, _vec, _norm, _lam) in toIter
        progress += 1
        @info @sprintf("%.3f", progress/length(toIter))
        hyperD = Dict{String,Real}("cancel"=>_cancel, "mean_k"=>_mean_k, "vec"=>_vec, "norm"=>_norm, "lam"=>_lam)
        network("WT.adj", "WP.adj"; header=true, o="net.bson", hyper=hyperD)
        isfile("sim_logFC.tsv") && rm("sim_logFC.tsv")
        logFC("net.bson", "sim_logFC.tsv")
        if isfile("sim_logFC.tsv")
            sim = loaddlm("sim_logFC.tsv"; header=true)[:, 2:end] |> Matrix
        else
            sim = zeros(size(ideal))
        end
        sim = vec(sim)
        hyperD["score"] = (sum(sim[ideal .== 1] .> thres) + sum(sim[ideal .== -1] .< -thres) + sum(abs.(sim[ideal .== 0]) .< thres)) / length(ideal)
        write(io, join([hyperD[k] for k in header], '\t') * '\n')
    end
end
