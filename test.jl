#!/usr/bin/env julia
cd("data")

include("PKTFX.jl")

PKTFX.network()
PKTFX.xgmml("WT.mat", "WP.mat", o="net.xgmml")

PKTFX.steadystate()
PKTFX.xgmml(o="steady.xgmml", X="steady_rpphi.mat")

PKTFX.iteratemodel(o="X_iter.mat")
PKTFX.xgmml("WT.mat", "WP.mat", o="iter.xgmml", X="X_iter.mat")

PKTFX.logFC("net.bson", o="X_sim.mat")
PKTFX.xgmml(o="sim.xgmml", X="X_sim.mat")

for i in 1:length(Inference.loss)
	println("conf=$i")
	# PKTFX.infer("X_iter.mat", 3, 3, "WT_iter_infer.mat", "WP_iter_infer.mat"; conf=i)
	# PKTFX.thres("WT_iter_infer.mat", "WT_iter_infer_thres.mat")
	# PKTFX.thres("WP_iter_infer.mat", "WP_iter_infer_thres.mat")
	# PKTFX.xgmml("WT_iter_infer_thres.mat", "WP_iter_infer_thres.mat", o="iter_infer_$i.xgmml")

	PKTFX.infer("X_sim.mat", 3, 3, "WT_sim_infer.mat", "WP_sim_infer.mat"; conf=i, lambda=1.)
	sleep(2)
	PKTFX.thres("WT_sim_infer.mat", "WT_sim_infer_thres.mat")
	PKTFX.thres("WP_sim_infer.mat", "WP_sim_infer_thres.mat")
	sleep(2)
	PKTFX.xgmml("WT_sim_infer_thres.mat", "WP_sim_infer_thres.mat", o="sim_infer_$i.xgmml")
end

# now open xgmml files in cytoscape and have a look
