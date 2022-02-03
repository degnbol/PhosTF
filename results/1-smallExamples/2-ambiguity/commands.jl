#!/usr/bin/env julia
@src "simulation/simulate";

nₜ = ncol(loaddlm("WT.adj"; header=true)) - 1
nₚ = ncol(loaddlm("WP.adj"; header=true)) - 1

network("WT.adj", "WP.adj"; header=true)
timeseries("net.bson")
plot_timeseries("sim.tsv", nₜ, nₚ; o="sim.pdf")
for mut_id in 1:nₜ+nₚ
    timeseries("net.bson", mut_id; u0="sim.tsv")
    plot_timeseries("sim_mut$mut_id.tsv", nₜ, nₚ; o="sim_mut$mut_id.pdf")
end
steadystate("net.bson")
for mut_id in 1:nₜ+nₚ
    steadystate("net.bson", mut_id; u0="sim.tsv")
end
logFC("net.bson", "sim_logFC.tsv")

