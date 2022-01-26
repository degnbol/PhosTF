#!/usr/bin/env julia
@src "simulation/simulate";

network("WT.adj", "WP.adj"; header=true)
timeseries("net.bson")
plot_timeseries("sim.tsv", 2, 3; o="sim.pdf")
for mut_id in 1:5
    #= timeseries("net.bson", mut_id; u0="sim.tsv") =#
    plot_timeseries("sim_mut$mut_id.tsv", 2, 3; o="sim_mut$mut_id.pdf")
end
steadystate("net.bson")
for mut_id in 1:5
    steadystate("net.bson", mut_id; u0="sim.tsv")
end
logFC("net.bson", "sim_logFC.tsv")


