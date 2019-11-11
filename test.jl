#!/usr/bin/env julia
cd("testdata/pres18c")
cd("../data_cas"); nₜ, nₚ = 3, 3
cd("../pres18c"); nₜ, nₚ = 3, 3
cd("../simi"); nₜ, nₚ = 2, 3
cd("../simi_neg"); nₜ, nₚ = 2, 3
cd("../cycl"); nₜ, nₚ = 2, 2
cd("../cycl_neg"); nₜ, nₚ = 2, 2


tryrm(fname) = try rm(fname) catch IOError end
# include("../../utilities/UnitTest.jl")
include("../../PKTFX.jl")

# PKTFX.network()
PKTFX.xgmml(o="net.xgmml")

PKTFX.simulate()
PKTFX.simulate(1; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")
PKTFX.simulate(2; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")
PKTFX.simulate(3; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")
PKTFX.simulate(4; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")

PKTFX.steadystate()
PKTFX.xgmml(X=["steady_r.mat", "steady_p.mat", "steady_phi.mat"], o="steady.xgmml")

tryrm("X_sim.mat"); PKTFX.logFC(o="X_sim.mat")
tryrm("sim.xgmml"); PKTFX.xgmml(X="X_sim.mat", o="sim.xgmml")

PKTFX.T("X_sim.mat", "T_sim.mat")
PKTFX.thres("T_sim.mat", "T_sim_thres.mat")
PKTFX.xgmml("T_sim_thres.mat", nₜ, nₚ; o="T.xgmml")

tryrm("WT_infer.mat"); tryrm("WP_infer.mat"); tryrm("WT_infer_thres.mat"); tryrm("WP_infer_thres.mat"); tryrm("sim_infer.xgmml")
PKTFX.infer("X_sim.mat", nₜ, nₚ; lambda=0., epochs=20000)
PKTFX.thres("WT_infer.mat", "WT_infer_thres.mat")
PKTFX.thres("WP_infer.mat", "WP_infer_thres.mat")
PKTFX.correct("WT_infer_thres.mat", "WP_infer_thres.mat"; ot="WT_infer_correct.mat", op="WP_infer_correct.mat")
PKTFX.xgmml("WT_infer_correct.mat", "WP_infer_correct.mat", o="sim_infer.xgmml")



# now open xgmml files in cytoscape and have a look
