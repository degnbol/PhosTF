#!/usr/bin/env julia

using Fire
include("PKTFX.jl")
# tryrm(fname) = try rm(fname) catch IOError end

@main function run(dataset::Integer)
    cd("$(dataset)_3")
    println(pwd())

    n, nₜ = size(PKTFX.loaddlm("sim/WT.mat"))
    _, nₚ = size(PKTFX.loaddlm("sim/WP.mat"))

    # tryrm("WT_infer.mat"); tryrm("WP_infer.mat"); tryrm("WT_infer_thres.mat"); tryrm("WP_infer_thres.mat")
    # tryrm("infer.xgmml")
    PKTFX.infer("sim/X_sim.mat", nₜ, nₚ; lambda=.1, epochs=5000)
    # PKTFX.infer("sim/X_sim.mat", nₜ, nₚ, "WT_infer_V.mat", "WP_infer_V.mat"; lambda=lambda, epochs=1000, PKPP="sim/PKPP.mat", WT="WT_infer.mat", WP="WP_infer.mat")
    # PKTFX.thres("WT_infer_V.mat", "WT_infer_thres.mat"; thres=1e-4)
    # PKTFX.thres("WP_infer_V.mat", "WP_infer_thres.mat"; thres=1e-4)
    # PKTFX.correct("WT_infer_thres.mat", "WP_infer_thres.mat")
end
