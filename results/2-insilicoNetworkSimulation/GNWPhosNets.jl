#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkSimulation
include(readchomp(`git root`) * "/src/simulation/simulate.jl");

mkpath("GNWPhosNets")

#= for n in [10, 100] =#
for n in [10]
    for i in 1:5
        for rep in 1:5
            println("net $n $i $rep")
            fnames = ["adjacencies/W$(s)_$(n)_$(i)-rep$(rep).adj" for s in "TP"]
            network(fnames...; header=true, o="GNWPhosNets/net_$(n)_$(i)-rep$(rep).bson")
        end
    end
end

