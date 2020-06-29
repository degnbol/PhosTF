#!/usr/bin/env julia

expanduser("~/cwd/simulation.jl") |> include

network()
simulate()
plot(3, 2; o="sim.pdf")
steadystate()
logFC()


