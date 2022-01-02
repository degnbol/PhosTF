#!/usr/bin/env julia
include(readchomp(`git root`) * "/src/simulation/simulate.jl");

network("WT.adj", "WP.adj"; header=true)
timeseries()
for mut_id in 1:4
    timeseries(mut_id; r0="sim_r.mat", p0="sim_p.mat", psi0="sim_psi.mat")
end
plot(2, 3; o="sim.pdf")
steadystate()
logFC()


