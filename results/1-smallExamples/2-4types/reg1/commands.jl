#!/usr/bin/env julia
# Run me from one of the folders that has the WT.adj and WP.adj
@src "simulation/simulate";

nₜ = ncol(loaddlm("WT.adj"; header=true)) - 1
nₚ = ncol(loaddlm("WP.adj"; header=true)) - 1

network("WT.adj", "WP.adj"; header=true, o="net.bson")

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

@src "inference/infer"
infer("sim_logFC.tsv", r"T.*", r"P.*"; epochs=10000, lambda_Bstar=0., lambda_absW=0.1)

@src "visualization/W2graph"
xgmml("net.bson"; o="net.xgmml")
xgmml("net.bson"; o="steady.xgmml", X="steady.tsv")
for mut_id in 1:nₜ+nₚ
    xgmml("net.bson"; o="steady_mut$mut_id.xgmml", X="steady_mut$mut_id.tsv", highlight=mut_id)
end
xgmml("net.bson"; o="ko.xgmml", X="sim_logFC.tsv")
xgmml("WT_infer.tsv", "WP_infer.tsv"; o="infer.xgmml")



