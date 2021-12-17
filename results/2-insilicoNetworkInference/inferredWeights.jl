#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkInference
ROOT=readchomp(`git root`)
include(ROOT * "/src/inference/infer.jl");
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet
using Glob

mkpath("inferredWeights")
mkpath("logs")

cd(ROOT * "/results/2-insilicoNetworkInference")
logFC_dir = glob("../*-insilicoNetworkSimulation/sim_logFCs/") |> only
#= logFC_fname = logFC_dir * "sim_logFC_100_1-rep1.mat" =#
WT = nothing
WP = nothing
WT_mask = nothing
WP_mask = nothing
TF = r"^TF[0-9]+$"
KP = r"^KP[0-9]+$"
mut_sep = nothing
#= infer(logFC_fname, TF, KP) =#

for n in [10]
    for (i, rep) in collect(product(1:5, 1:5))
        SUF="_$(n)_$(i)-rep$(rep)"
        logFC_fname = logFC_dir * "sim_logFC$SUF.mat"
        infer(logFC_fname, TF, KP, "inferredWeights/WT$SUF.adj", "inferredWeights/WP$SUF.adj"; log="logs/log$SUF.tsv")
    end
end

