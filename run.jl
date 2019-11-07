#!/usr/bin/env julia
include("PKTFX.jl")
tryrm(fname) = try rm(fname) catch IOError end

nₜ, nₚ = 3, 3

for i in 1:2length(Inference.loss) tryrm("sim_infer_$i.xgmml") end
for i in 1:length(Inference.loss)
	println("conf=$i")
	tryrm("WT_sim_infer.mat"); tryrm("WP_sim_infer.mat")
	PKTFX.infer("X_sim.mat", nₜ, nₚ, "WT_sim_infer.mat", "WP_sim_infer.mat"; conf=i, lambda=.1, epochs=4000)
	sleep(2)
	PKTFX.thres("WT_sim_infer.mat", "WT_sim_infer_thres.mat")
	PKTFX.thres("WP_sim_infer.mat", "WP_sim_infer_thres.mat")
	sleep(2)
	PKTFX.xgmml("WT_sim_infer_thres.mat", "WP_sim_infer_thres.mat", o="sim_infer_$i.xgmml")
	tryrm("WT_sim_infer_thres.mat"); tryrm("WP_sim_infer_thres.mat")
end
