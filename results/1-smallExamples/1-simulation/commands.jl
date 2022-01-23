#!/usr/bin/env julia
@src "simulation/simulate";

network("WT.adj", "WP.adj"; header=true)
timeseries("net.bson")
for mut_id in 1:5
    timeseries("net.bson", mut_id; t0="sim.tsv")
end
plot("sim.tsv", 2, 3; o="sim.pdf")
steadystate("net.bson")
logFC("net.bson")


