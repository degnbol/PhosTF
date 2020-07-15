#!/usr/bin/env julia

expanduser("~/cwd/simulation.jl") |> include

network()
simulate()
for mut_id in 1:4
    simulate(mut_id; r0="sim_r.mat", p0="sim_p.mat", psi0="sim_psi.mat")
end
plot(3, 2; o="sim.pdf")
steadystate()
logFC()


