#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkSimulation
include(readchomp(`git root`) * "/src/simulation/simulate.jl");

mkpath("GNWPhosNets")

for n in [10, 100]
    for i in 1:5
        for rep in 1:5
            fnames = ["adjacencies/W$(s)_$(n)_$(i)-rep$(rep).adj" for s in "TP"]
            network(fnames...; header=true, o="GNWPhosNets/net_$(n)_$(i)-rep$(rep).bson")
        end
    end
end

