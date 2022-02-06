#!/usr/bin/env julia
include(readchomp(`git root`) * "/src/simulation/simulate.jl");
include(readchomp(`git root`) * "/src/inference/infer.jl");
using .Threads: @threads
using .Iterators: product # since @threads does not support nested loops yet
using Glob

nzfile(fname) = isfile(fname) && filesize(fname) > 0

mkpath("GNWPhosNets")
mkpath("sim_logFCs")
mkpath("inferredWeights")
mkpath("logs")
mkpath("aucs")

n = 10
rep = 1
@threads for (i, _cancel, _mean_k, _vec, _norm, _lam) in collect(product(1:5, [false, true], [true, false], [true, false], 0:2, 0:4))
    hyper = "$_cancel-$_mean_k-$_vec-$_norm-$_lam"
    hyperD = Dict("cancel"=>_cancel, "mean_k"=>_mean_k, "vec"=>_vec, "norm"=>_norm, "lam"=>_lam)
    SUF="_$(n)_$(i)-rep$(rep)"
    println("net $n $i $rep $hyper")
    netdir = "GNWPhosNets/$hyper" 
    mkpath(netdir)
    fnames = ["../2-insilicoNetworkSimulation/adjacencies/W$(s)$SUF.adj" for s in "TP"]
        nzfile("$netdir/net$SUF.bson") || network(fnames...; header=true, o="$netdir/net$SUF.bson", hyper=hyperD)

    println("sim $n $i $rep $hyper")
    simdir = "sim_logFCs/$hyper" 
    mkpath(simdir)
    nzfile("$simdir/sim_logFC$SUF.tsv") || logFC("$netdir/net$SUF.bson", "$simdir/sim_logFC$SUF.tsv");

    println("inf $n $i $rep $hyper")
    logdir = "logs/$hyper"
    infdir = "inferredWeights/$hyper"
    mkpath(logdir)
    mkpath(infdir)
    for λBstar ∈ [0, .1, .5, 1.]
        for λabsW ∈ [0, .1, .5, 1.]
            ot = "$infdir/WT$SUF-$λBstar-$λabsW.adj"
            op = "$infdir/WP$SUF-$λBstar-$λabsW.adj"
            (nzfile(ot) && nzfile(op)) || infer("$simdir/sim_logFC$SUF.tsv", r"^T[0-9]+$", r"^P[0-9]+$", "$infdir/WT$SUF-$λBstar-$λabsW.adj", "$infdir/WP$SUF-$λBstar-$λabsW.adj"; log="$logdir/log$SUF-$λBstar-$λabsW.tsv", lambda_Bstar=λBstar, lambda_absW=λabsW)
        end
    end
end

