#!/usr/bin/env julia
include("PKTFX.jl")
tryrm(fname) = try rm(fname) catch IOError end

# run from a folder with a name like "lambda_01_001_050"

lambda_W, lambda_B = split(basename(pwd()), '_')[2:3]
lambda_W = parse(Float64, '.' * lambda_W[2:end])
lambda_B = parse(Float64, '.' * lambda_B[2:end])

for net ∈ 1:5
    cd("$(net)_1")
    println(pwd())

    n, nₜ = size(PKTFX.loaddlm("sim/WT.mat"))
    _, nₚ = size(PKTFX.loaddlm("sim/WP.mat"))

    PKTFX.infer("sim/X_sim.mat", nₜ, nₚ; lambda_W=lambda_W, lambda_B=lambda_B, epochs=10000)
    PKTFX.infer("sim/X_sim.mat", nₜ, nₚ, "WT_infer_V.mat", "WP_infer_V.mat"; lambda_W=lambda_W, lambda_B=lambda_B, epochs=10000, PKPP="sim/PKPP.mat", WT="WT_infer.mat", WP="WP_infer.mat")
    PKTFX.thres("WT_infer_V.mat", "WT_infer_thres.mat"; thres=1e-4)
    PKTFX.thres("WP_infer_V.mat", "WP_infer_thres.mat"; thres=1e-4)
    PKTFX.correct("WT_infer_thres.mat", "WP_infer_thres.mat")

    cd("..")
end


