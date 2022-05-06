#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkInference
@src "inference/infer"
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet
using Glob

mkpath("inferredWeights")
mkpath("logs")

WT = nothing
WP = nothing
WT_mask = nothing
WP_mask = nothing
TF = r"^TF?[0-9]+$"
KP = r"^K?P[0-9]+$"
mut_sep = nothing

logFC_dir = glob("../*-insilicoNetworkSimulation/sim_logFCs/") |> only
# logFC_fname = logFC_dir * "sim_logFC_100_1-rep1.tsv"
# infer(logFC_fname, TF, KP)

for n in [10, 100]
    @threads for (i, rep) in collect(product(1:5, 1:5))
        SUF="_$(n)_$(i)-rep$(rep)"
        logFC_fname = logFC_dir * "sim_logFC$SUF.tsv"
        log_fname = "logs/log$SUF.tsv"
        infer(logFC_fname, TF, KP, "inferredWeights/WT$SUF.adj", "inferredWeights/WP$SUF.adj"; log=log_fname, lambda_Bstar=0.1, lambda_absW=0., epochs=20000)
    end
end

