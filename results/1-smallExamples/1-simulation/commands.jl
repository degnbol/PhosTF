#!/usr/bin/env julia
include(readchomp(`git root`) * "/src/simulation/simulate.jl");

network("WT.adj", "WP.adj"; header=true)
timeseries("net.bson")
for mut_id in 1:5
    timeseries("net.bson", mut_id; t0="sim.tsv")
end
plot(2, 3; o="sim.pdf")
steadystate("net.bson")
logFC("net.bson")


