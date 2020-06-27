#!/usr/bin/env julia
expanduser("~/cwd/inference.jl") |> include;

for dir in ARGS
    println(dir); cd(dir)

    nᵥ, nₜ = size(loaddlm("WT.mat"))
    _,  nₚ = size(loaddlm("WP.mat"))

    infer("X_sim.mat", nₜ, nₚ; epochs=15000)
        
    cd("..")
end

