#!/usr/bin/env julia

expanduser("~/cwd/cytoscape.jl") |> include

xgmml(; o="net.xgmml")
xgmml(; o="steady.xgmml", X=["steady_r.mat", "steady_p.mat", "steady_psi.mat"])
xgmml(; o="sim.xgmml", X="X_sim.mat")
xgmml("WT_infer.mat", "WP_infer.mat"; o="infer.xgmml")


