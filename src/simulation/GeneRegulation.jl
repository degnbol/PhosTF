#!/usr/bin/env julia
isdefined(Main, :Model) || include("../inference/Model.jl")

"""
Structs with data defining a gene regulation network and its regulation mechanisms. All functions for initializing the values are found here, including randomized initialization.
"""
module GeneRegulation
using Statistics: mean
import JSON3
import ..Model: nₒnₜnₚ, WₜWₚ

export Network, Gene
export drdt, dpdt, dψdt
export nₒnₜnₚ, estimate_Wₜ


const weak_activation = .25
const strong_activation = .9
const noise_activation = .25weak_activation


include("RegulatoryModule.jl")
include("Gene.jl")
include("phos_edges.jl")
include("Network.jl")

# define struct types for JSON3 to be able to read/write them
JSON3.StructType(::Type{RegulatoryModule}) = JSON3.Struct()
JSON3.StructType(::Type{Gene}) = JSON3.Struct()
JSON3.StructType(::Type{Network}) = JSON3.Struct()


nₒnₜnₚ(net::Network) = net.nₒ,net.nₜ,net.nₚ


random_t½() = TruncNormal(5, 50)
"""
We have exponential decay, the half-life and the decay rate are thus related by:
t½ = ln(2) / λ ⟹
λ = ln(2) / t½
"""
random_λ(n::Int) = log(2) ./ rand(random_t½(), n)
"""
For λ₊, λ₋ to have lower values for nodes that are mostly negatively regulated.
"""
function random_λ(mat::Matrix)
    # weigh by the fraction of regulators that regulate positively.
    positives = sum(mat .> 0; dims=2) |> vec
    negatives = sum(mat .< 0; dims=2) |> vec
    λ₊ = random_λ(size(mat,1)) .* negatives ./ (positives .+ negatives)
    λ₋ = random_λ(size(mat,1)) .* positives ./ (positives .+ negatives)
    # NaN from div zero which means there are no regulators of a node. In that case it's activation is static.
    noreg = positives .+ negatives .== 0
    λ₊[noreg] .= random_λ(sum(noreg))
    λ₋[noreg] .= 0
    λ₊, λ₋
end


"""
Mean activation μ given ψ.
ψ: 1D array. Active nondim concentrations.
"""
function μ(m::RegulatoryModule, ψ::AbstractVector{<:AbstractFloat})
	χ = (ψ[m.inputs] ./ m.k) .^ m.ν
	activator_prod = prod(χ[1:m.n_activators])
	if m.complex
		denom = 1 + activator_prod
		if m.n_repressors > 0 denom += prod(χ) end
	else
		denom = prod(1 .+ χ)
	end
	activator_prod / denom
end
μ(ms::Vector{RegulatoryModule}, ψ::AbstractVector{<:AbstractFloat}) = [μ(m, ψ) for m in ms]
"""
Fraction of max activation for a given gene when active TFs are found at a given concentration.
"""
function f(gene::Gene, ψ::AbstractVector{<:AbstractFloat})
    isempty(gene.modules) && 1
	μs = μ(gene.modules, ψ)
	# get P{state} for all states, each state is a unique combination of modules
	P = [prod(μs[state]) * prod(1 .- μs[.!state]) for state in states(length(gene.modules))]
	sum(gene.α .* P)
end
f(genes::Vector{Gene}, ψ::AbstractVector{<:AbstractFloat}) = [f(gene, ψ) for gene in genes]


"Estimate the effect on f for all genes when a given TF has either ψ=weak or ψ=strong."
function estimate_Wₜ(net::Network, i::Integer, basal_activation::AbstractFloat)
	basal = fill(basal_activation, net.nᵥ)
	ψ = copy(basal); ψ[i] = 1
	f(net.genes, ψ) - f(net.genes, basal)
end
function estimate_Wₜ(net::Network, i::Integer)
	mean(estimate_Wₜ(net, i, activation)
	for activation in [0, noise_activation, weak_activation])
end
"Estimate Wₜ by comparing the effect on f when any TF has ψ=weak or ψ=strong."
estimate_Wₜ(net::Network) = hcat([estimate_Wₜ(net, i) for i in net.nₚ+1:net.nₚ+net.nₜ]...)

end;

