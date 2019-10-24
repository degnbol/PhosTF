#!/usr/bin/env julia
cd("testdata/data_cas")

include("PKTFX.jl")

PKTFX.network()
PKTFX.xgmml(o="net.xgmml")

PKTFX.simulate()
PKTFX.plot(o="tmp.pdf")

PKTFX.steadystate()
run(`cp steady_phi.mat steady_phi0.mat`)  # copy phi values to file where we will append a single zero
run(pipeline(`echo 0.0`, stdout="steady_phi0.mat", append=true))  # append a single zero bc nₓ==1
run(pipeline(`paste -d" " steady_r.mat steady_p.mat steady_phi0.mat`, "steady_rpphi.mat"))
PKTFX.xgmml(X="steady_rpphi.mat", o="steady.xgmml")

# PKTFX.iteratemodel(o="X_iter.mat")
# PKTFX.xgmml("WT.mat", "WP.mat", o="iter.xgmml", X="X_iter.mat")
rm("X_sim.mat"); PKTFX.logFC(o="X_sim.mat")
rm("sim.xgmml"); PKTFX.xgmml(X="X_sim.mat", o="sim.xgmml")

for i in 1:length(Inference.loss)
	println("conf=$i")
	# PKTFX.infer("X_iter.mat", 3, 3, "WT_iter_infer.mat", "WP_iter_infer.mat"; conf=i)
	# PKTFX.thres("WT_iter_infer.mat", "WT_iter_infer_thres.mat")
	# PKTFX.thres("WP_iter_infer.mat", "WP_iter_infer_thres.mat")
	# PKTFX.xgmml("WT_iter_infer_thres.mat", "WP_iter_infer_thres.mat", o="iter_infer_$i.xgmml")

	PKTFX.infer("X_sim.mat", 3, 3, "WT_sim_infer.mat", "WP_sim_infer.mat"; conf=i, lambda=.1)
	sleep(2)
	PKTFX.thres("WT_sim_infer.mat", "WT_sim_infer_thres.mat")
	PKTFX.thres("WP_sim_infer.mat", "WP_sim_infer_thres.mat")
	sleep(2)
	PKTFX.xgmml("WT_sim_infer_thres.mat", "WP_sim_infer_thres.mat", o="sim_infer_$i.xgmml")
end

# now open xgmml files in cytoscape and have a look


