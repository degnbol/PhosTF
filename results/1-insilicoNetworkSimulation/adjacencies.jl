#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkSimulation
ROOT=readchomp(`git root`)
include(ROOT * "/src/simulation/weight.jl");

mkpath("adjacencies")

for n in [10, 100]
    for i in 1:5
        for rep in 1:5
            fname = ROOT * "/data/DREAM4/adjacencies/goldstandard_$(n)_$(i).adj"
            random(fname, Int(.3 * n); WT="WT_$(n)_$(i)-rep$(rep).adj", WP="WP_$(n)_$(i)-rep$(rep).adj", pretty=true);
        end
    end
end
