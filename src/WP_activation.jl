#!/usr/bin/env julia
include("simulation.jl")

for net ∈ 1:5 for sample ∈ 1:5
    cd("$(net)_$(sample)")
    WP_activation(o="WP_activation.mat")
    cd("..")
end end
