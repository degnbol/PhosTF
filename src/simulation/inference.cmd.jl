#!/usr/bin/env julia
include("../../inference.jl");

"""
Inference runs on multiple simulation folders given in ARGS.
"""

for dir in ARGS
    println(dir); cd(dir)

    nᵥ, nₜ = size(loaddlm("WT.mat"))
    _,  nₚ = size(loaddlm("WP.mat"))

    open("infer.log", "w") do log redirect_stdout(log) do
            infer("X_sim.mat", nₜ, nₚ; epochs=15000, lambda=1.)
    end end
    
    rm("WT.tmp.mat")
    rm("WP.tmp.mat")

    cd("..")
end

