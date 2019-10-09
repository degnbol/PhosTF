#!/usr/bin/env julia
include("Model.jl")
include("GeneRegulationNetwork.jl")

"""
Module intended to hold all structs with data defining our gene regulation network and its regulation mechanisms. All functions for initializing the values are found here, including randomized initialization.
"""
module GeneRegulation
using Statistics: mean
import JSON3
import ..Model: nₓnₜnₚ
using ..GeneRegulationNetwork
using ..GeneRegulationNetwork.GeneRegulationGene

export Network, Gene
export drdt, dpdt, dϕdt
export nₓnₜnₚ, estimate_Wₜ

# define struct types for JSON3 to be able to read/write them
JSON3.StructType(::Type{GeneRegulationGene.RegulatoryModule}) = JSON3.Struct()
JSON3.StructType(::Type{Gene}) = JSON3.Struct()
JSON3.StructType(::Type{Network}) = JSON3.Struct()


"Estimate the effect on f for all genes when a given TF has either ϕ=weak or ϕ=strong."
function estimate_Wₜ(net::Network, i::Integer, basal_activation::AbstractFloat)
	basal = fill(basal_activation, net.n)
	ψ = copy(basal); ψ[i] = 1
	f(net.genes, ψ) - f(net.genes, basal)
end
function estimate_Wₜ(net::Network, i::Integer)
	mean(estimate_Wₜ(net, i, activation)
	for activation in [0, GeneRegulationGene.noise_activation, GeneRegulationGene.weak_activation])
end
"Estimate Wₜ by comparing the effect on f when any TF has ϕ=weak or ϕ=strong."
estimate_Wₜ(net::Network) = hcat([estimate_Wₜ(net, i) for i in net.nₚ+1:net.nₚ+net.nₜ]...)


function nₓnₜnₚ(net::Network)
	println("yay")
	net.n-(net.nₜ+net.nₚ),net.nₜ,net.nₚ
end

end;

