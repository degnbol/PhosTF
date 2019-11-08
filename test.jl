#!/usr/bin/env julia
cd("testdata/pres18c")
cd("../data_cas")
cd("../pres18c")
nₜ, nₚ = 3, 3


tryrm(fname) = try rm(fname) catch IOError end
# include("../../utilities/UnitTest.jl")
include("../../PKTFX.jl")

PKTFX.network()
PKTFX.xgmml(o="net.xgmml")

PKTFX.simulate()
PKTFX.simulate(2; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")
PKTFX.simulate(4; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")

PKTFX.steadystate()
PKTFX.xgmml(X=["steady_r.mat", "steady_p.mat", "steady_phi.mat"], o="steady.xgmml")

tryrm("X_sim.mat"); PKTFX.logFC(o="X_sim.mat")
tryrm("sim.xgmml"); PKTFX.xgmml(X="X_sim.mat", o="sim.xgmml")

tryrm("sim_infer.xgmml")
tryrm("WT_sim_infer.mat"); tryrm("WP_sim_infer.mat")
PKTFX.infer("X_sim.mat", nₜ, nₚ, "WT_sim_infer.mat", "WP_sim_infer.mat"; lambda=10, epochs=10000)
sleep(1)
PKTFX.thres("WT_sim_infer.mat", "WT_sim_infer_thres.mat")
PKTFX.thres("WP_sim_infer.mat", "WP_sim_infer_thres.mat")
sleep(1)
PKTFX.xgmml("WT_sim_infer_thres.mat", "WP_sim_infer_thres.mat", o="sim_infer.xgmml")
tryrm("WT_sim_infer_thres.mat"); tryrm("WP_sim_infer_thres.mat")


# now open xgmml files in cytoscape and have a look
