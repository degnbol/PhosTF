#!/usr/bin/env julia
# cd("testdata/pres18c")
# cd("../data_cas"); nₜ, nₚ = 3, 3
# cd("../pres18c"); nₜ, nₚ = 3, 3
# cd("../simi"); nₜ, nₚ = 2, 3
# cd("../simi_neg"); nₜ, nₚ = 2, 3
# cd("../cycl"); nₜ, nₚ = 2, 2
# cd("../cycl_neg"); nₜ, nₚ = 2, 2
using Test
include("PKTFX.jl")
tryrm(fname) = try rm(fname) catch IOError end
# include("../../utilities/UnitTest.jl")


# for attempt ∈ 1:10
#     println("attempt: $attempt")
#     PKTFX.network()
#     tryrm("X_sim.mat")
#     try
#         # test that there are no warnings
#         @test_logs (:info, "logFC values simulated") PKTFX.logFC(o="X_sim.mat")
#     catch
#         continue  # attempt again
#     end
#     break
# end

# PKTFX.xgmml(o="net.xgmml")

# PKTFX.simulate()
# PKTFX.simulate(1; r0="sim_r.mat", p0="sim_p.mat", phi0="sim_phi.mat")

# PKTFX.steadystate()
# PKTFX.xgmml(X=["steady_r.mat", "steady_p.mat", "steady_phi.mat"], o="steady.xgmml")

n, nₜ = size(PKTFX.loaddlm("WT.mat"))
_, nₚ = size(PKTFX.loaddlm("WP.mat"))

# tryrm("sim.xgmml"); PKTFX.xgmml(X="X_sim.mat", o="sim.xgmml")

# PKTFX.T("X_sim.mat", "T_sim.mat")
# PKTFX.thres("T_sim.mat", "T_sim_thres.mat")
# PKTFX.xgmml("T_sim_thres.mat", nₜ, nₚ; o="T.xgmml")

tryrm("WT_infer.mat"); tryrm("WP_infer.mat"); tryrm("WT_infer_thres.mat")
tryrm("WP_infer_thres.mat"); tryrm("infer.xgmml")
PKTFX.infer("X_sim.mat", nₜ, nₚ; lambda_W=0.01, lambda_B=0.1, epochs=10000)
PKTFX.infer("X_sim.mat", nₜ, nₚ; lambda_W=0., lambda_B=1., epochs=10000, PKPP="PKPP.mat", WT="WT_infer.mat", WP="WP_infer.mat")
PKTFX.thres("WT_infer.mat", "WT_infer_thres.mat")
PKTFX.thres("WP_infer.mat", "WP_infer_thres.mat")
PKTFX.xgmml("WT_infer_thres.mat", "WP_infer_thres.mat", o="infer.xgmml")

# now open xgmml files in cytoscape and have a look

