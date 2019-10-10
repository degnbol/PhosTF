#!/usr/bin/env julia
include("../Model.jl")

"""
Structs with data defining a gene regulation network and its regulation mechanisms. All functions for initializing the values are found here, including randomized initialization.
"""
module GeneRegulation
using Statistics: mean
import JSON3
import ..Model: nₓnₜnₚ

export Network, Gene
export drdt, dpdt, dϕdt
export nₓnₜnₚ, estimate_Wₜ


const weak_activation = .25
const noise_activation = .25weak_activation


include("regulatory_module.jl")
include("gene.jl")
include("network.jl")

# define struct types for JSON3 to be able to read/write them
JSON3.StructType(::Type{RegulatoryModule}) = JSON3.Struct()
JSON3.StructType(::Type{Gene}) = JSON3.Struct()
JSON3.StructType(::Type{Network}) = JSON3.Struct()


nₓnₜnₚ(net::Network) = net.nₓ,net.nₜ,net.nₚ


"""
Concentration of active protein, which is either phosphorylated or unphosphorylated concentration of the protein, depending on a bool.
"""
ψ(p, ϕ, phos_activation) = @. phos_activation * ϕ + (!phos_activation) * (p-ϕ)

"""
Mean activation μ given ψ.
ψ: 1D array. Active nondim concentrations.
"""
function μ(m::RegulatoryModule, ψ::Vector{<:AbstractFloat})
	χ = (ψ[m.inputs] ./ m.k) .^ m.ν
	activator_prod = prod(χ[1:m.n_activators])
	if m.complex
		denom = 1 .+ activator_prod
		if m.n_repressors > 0 denom += prod(χ) end
	else
		denom = prod(1 .+ χ)
	end
	activator_prod / denom
end
μ(ms::Vector{RegulatoryModule}, ψ::Vector{<:AbstractFloat}) = [μ(m, ψ) for m in ms]
"""
Fraction of max activation for a given gene when active TFs are found at a given concentration.
"""
function f(gene::Gene, ψ::Vector{<:AbstractFloat})
	if isempty(gene.modules) return 1 end
	μs = μ(gene.modules, ψ)
	# get P{state} for all states, each state is a unique combination of modules
	P = [prod(μs[state]) * prod(1 .- μs[.!state]) for state in states(length(gene.modules))]
	sum(gene.α .* P)
end
f(genes::Vector{Gene}, ψ::Vector{<:AbstractFloat}) = [f(gene, ψ) for gene in genes]


"Estimate the effect on f for all genes when a given TF has either ϕ=weak or ϕ=strong."
function estimate_Wₜ(net::Network, i::Integer, basal_activation::AbstractFloat)
	basal = fill(basal_activation, net.n)
	ψ = copy(basal); ψ[i] = 1
	f(net.genes, ψ) - f(net.genes, basal)
end
function estimate_Wₜ(net::Network, i::Integer)
	mean(estimate_Wₜ(net, i, activation)
	for activation in [0, noise_activation, weak_activation])
end
"Estimate Wₜ by comparing the effect on f when any TF has ϕ=weak or ϕ=strong."
estimate_Wₜ(net::Network) = hcat([estimate_Wₜ(net, i) for i in net.nₚ+1:net.nₚ+net.nₜ]...)


end;

