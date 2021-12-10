#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkConstruction
ROOT=readchomp(`git root`)
# comment this out to use your own env:
using Pkg; Pkg.activate(ROOT * "/src")
include(ROOT * "/src/simulation/simulate.jl");
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet

mkpath("sim_logFCs")

for n in [10, 100]
    @threads for (i, rep) in collect(product(1:5, 1:5))
        logFC("GNWPhosNets/net_$(n)_$(i)-rep$(rep).bson"; o="sim_logFCs/sim_logFC_$(n)_$(i)-rep$(rep).mat");
    end
end

