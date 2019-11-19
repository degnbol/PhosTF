#!/usr/bin/env julia
include("PKTFX.jl")

for net ∈ 1:5 for sample ∈ 1:5
    cd("$(net)_$(sample)")
    PKTFX.PKPP("WP.mat", "PKPP.mat")
    cd("..")
end end
