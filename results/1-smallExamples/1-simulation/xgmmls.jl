#!/usr/bin/env julia
include(readchomp(`git root`) * "/src/visualization/W2graph.jl");

xgmml("net.bson"; o="net.xgmml")
xgmml("net.bson"; o="steady.xgmml", X=["steady_r.mat", "steady_p.mat", "steady_psi.mat"])
xgmml("net.bson"; o="ko.xgmml", X="X_sim.mat")
xgmml("WT_infer.mat", "WP_infer.mat"; o="infer.xgmml")


