#!/usr/bin/env julia
"""
Structs with data defining a gene regulation network and its regulation mechanisms.
All functions for initializing the values are found here, including randomized initialization.
"""
module GeneRegulation
using Statistics: mean
import JSON3

export Network, Gene
export drdt, dpdt, dψdt
export estimate_Wₜ

include("Network.jl")

# define struct types for JSON3 to be able to read/write them
JSON3.StructType(::Type{RegulatoryModule}) = JSON3.Struct()
JSON3.StructType(::Type{Gene}) = JSON3.Struct()
JSON3.StructType(::Type{Network}) = JSON3.Struct()


"Estimate the effect on f for all genes when a given TF has either ψ=weak or ψ=strong."
function estimate_Wₜ(net::Network, i::Integer, basal_activation::AbstractFloat)
	basal = fill(basal_activation, net.nᵥ)
	ψ = copy(basal); ψ[i] = 1
    f.(net.genes, Ref(ψ)) .- f.(net.genes, Ref(basal))
end
function estimate_Wₜ(net::Network, i::Integer)
    noise_activation = .25weak_activation
	mean(estimate_Wₜ(net, i, activation) for activation in [0, noise_activation, weak_activation])
end
"""
Estimate Wₜ by comparing the effect on f when any TF has ψ=weak or ψ=strong.
Assumes nodes are ordered with TFs first.
"""
estimate_Wₜ(net::Network) = hcat([estimate_Wₜ(net, i) for i in 1:net.nₜ]...)

end;

