#!/usr/bin/env julia
cd("testdata/data_cas")
tryrm(fname) = try rm(fname) catch IOError end
include("../../utilities/UnitTest.jl")
include("../../PKTFX.jl")

nₜ, nₚ = 3, 3

# PKTFX.network()
# PKTFX.xgmml(o="net.xgmml")

# PKTFX.simulate()
# PKTFX.plot(o="tmp.pdf")

# PKTFX.steadystate()
# PKTFX.xgmml(X=["steady_r.mat", "steady_p.mat", "steady_phi.mat"], o="steady.xgmml")

# tryrm("X_sim.mat"); PKTFX.logFC(o="X_sim.mat")
# tryrm("sim.xgmml"); PKTFX.xgmml(X="X_sim.mat", o="sim.xgmml")
# PKTFX.swapPT("X_sim.mat", "X_sim_TP.mat"; n1=3, n2=3) # in case you wanna use the python code

for i in 1:2length(Inference.loss) tryrm("sim_infer_$i.xgmml") end
for i in 1:length(Inference.loss)
	println("conf=$i")
	tryrm("WT_sim_infer.mat"); tryrm("WP_sim_infer.mat")
	PKTFX.infer("X_sim.mat", nₜ, nₚ, "WT_sim_infer.mat", "WP_sim_infer.mat"; conf=i, lambda=.1)
	sleep(2)
	PKTFX.thres("WT_sim_infer.mat", "WT_sim_infer_thres.mat")
	PKTFX.thres("WP_sim_infer.mat", "WP_sim_infer_thres.mat")
	sleep(2)
	PKTFX.xgmml("WT_sim_infer_thres.mat", "WP_sim_infer_thres.mat", o="sim_infer_$i.xgmml")
	tryrm("WT_sim_infer_thres.mat"); tryrm("WP_sim_infer_thres.mat")
end

begin # an attempt at using the B method 
	for i in 1:2length(Inference.loss_B) tryrm("sim_inferB_$i.xgmml") end
	for i in 1:length(Inference.loss_B)
		println("conf=$i")
		tryrm("WT_sim_inferB.mat"); tryrm("WP_sim_inferB.mat")
		PKTFX.inferB("B_sim.mat", nₚ, "WT_sim_inferB.mat", "WP_sim_inferB.mat"; conf=i, lambda=.1)
		sleep(2)
		PKTFX.thres("WT_sim_inferB.mat", "WT_sim_inferB_thres.mat")
		PKTFX.thres("WP_sim_inferB.mat", "WP_sim_inferB_thres.mat")
		sleep(2)
		PKTFX.xgmml("WT_sim_inferB_thres.mat", "WP_sim_inferB_thres.mat", o="sim_inferB_$i.xgmml")
		tryrm("WT_sim_inferB_thres.mat"); tryrm("WP_sim_inferB_thres.mat")
	end
end

# now open xgmml files in cytoscape and have a look

begin # for using the old python code
	PKTFX.swapPT("W_TP.mat", "W_py.mat"; n1=3, n2=3)
	PKTFX.thres("W_py.mat", "W_py_thres.mat")
	wt, wp = Model.WₜWₚ(PKTFX.loaddlm("W_py_thres.mat"), 3, 3)
	PKTFX.savedlm("WT_py_thres.mat", wt)
	PKTFX.savedlm("WP_py_thres.mat", wp)
	PKTFX.xgmml("WT_py_thres.mat", "WP_py_thres.mat", o="py_infer.xgmml")
end

