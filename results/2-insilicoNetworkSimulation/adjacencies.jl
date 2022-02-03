#!/usr/bin/env julia
# Run from `git root`/results/*-insilicoNetworkSimulation
ROOT=readchomp(`git root`)
include("$ROOT/src/simulation/weight.jl");

mkpath("adjacencies")

for n in [10, 100]
    for i in 1:5
        for rep in 1:5
            println("$n $i $rep")
            infname = ROOT * "/data/DREAM4/adjacencies/goldstandard_$(n)_$(i).adj"
            outfnames = ["adjacencies/W$(s)_$(n)_$(i)-rep$(rep).adj" for s in "TP"]
            random(infname, Int(.3 * n); WT=outfnames[1], WP=outfnames[2], pretty=true);
            validate(outfnames...)
        end
    end
end
