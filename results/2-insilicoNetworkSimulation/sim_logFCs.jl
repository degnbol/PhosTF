#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkSimulation
include(readchomp(`git root`) * "/src/simulation/simulate.jl");
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet

mkpath("sim_logFCs")

println("beginning sim")
#= for n in [10, 100] =#
for n in [10]
    @threads for (i, rep) in collect(product(1:5, 1:5))
        println("sim $n $i $rep")
        logFC("GNWPhosNets/net_$(n)_$(i)-rep$(rep).bson", "sim_logFCs/sim_logFC_$(n)_$(i)-rep$(rep).tsv");
    end
end

