#!/usr/bin/env julia
ROOT=readchomp(`git root`)
# comment this out to use your own env:
using Pkg; Pkg.activate(ROOT * "/src")
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
