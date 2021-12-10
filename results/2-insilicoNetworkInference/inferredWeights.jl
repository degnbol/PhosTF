#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkInference
ROOT=readchomp(`git root`)
include(ROOT * "/src/inference/infer.jl");
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet
using Glob

mkpath("inferredWeights")

cd(ROOT * "/results/2-insilicoNetworkInference")
logFC_dir = glob("../*-insilicoNetworkSimulation/sim_logFCs/") |> only
X = logFC_dir * "sim_logFC_10_1-rep1.mat"
J = nothing
WT = nothing
WP = nothing
WT_mask = nothing
WP_mask = nothing
nₜ = 7
nₚ = 3

infer(X, 7, 3)

#= for n in [10, 100] =#
#=     @threads for (i, rep) in collect(product(1:5, 1:5)) =#
#=          =#
#=     end =#
#= end =#

