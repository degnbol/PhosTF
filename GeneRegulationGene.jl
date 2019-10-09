#!/usr/bin/env julia
include("utilities/ArrayUtils.jl")

"""
Gene and regulatory module structs for gene regulation simulation.
"""
module GeneRegulationGene
using Distributions: Uniform, TruncatedNormal
using ..ArrayUtils

export Gene
export f, ψ

const weak_activation = .25
const noise_activation = .25weak_activation

"""
State s is interpreted as a binary number, where bit k indicates whether module k is active (True) or inactive (False) in this state.
n: number of modules to combine in each unique way
"""
states(n) = (binary(i-1, n) for i in 1:2^n)


struct RegulatoryModule
	n_activators::Int
	n_repressors::Int
	inputs::Vector{Int}  # int indexes for activators then repressors among all proteins
	ν::Vector{Float64}
	k::Vector{Float64}
	complex::Bool
	inhibitor::Bool
	# constructor for JSON3
	function RegulatoryModule(n_activators::Integer, n_repressors::Integer, inputs::Vector{<:Integer}, ν::Vector{<:AbstractFloat}, k::Vector{<:AbstractFloat}, complex::Bool, inhibitor::Bool)
		new(n_activators, n_repressors, inputs, ν, k, complex, inhibitor)
	end
	function RegulatoryModule(activators::Vector, repressors::Vector, ν::Vector, k::Vector, complex::Bool, inhibitor::Bool)
		new(length(activators), length(repressors), vcat(activators, repressors), ν, k, complex, inhibitor)
	end
	function RegulatoryModule(activators::Vector, repressors::Vector)
		n = length(activators) + length(repressors)
		ν = random_ν(n)
		k = random_k(n)
		complex = rand([true, false])
		inhibitor = random_inhibitor(length(activators), length(repressors))
		if inhibitor
			activators, repressors = repressors, activators
		end
		RegulatoryModule(activators, repressors, ν, k, complex, inhibitor)
	end
	
	"""
	Hill coeficient ν for each input protein.
	"""
	random_ν(n::Integer) = rand(TruncatedNormal(2, 2, 1, 10), n)
	"""
	Dissociation constant k for each input protein.
	"""
	random_k(n::Integer) = rand(Uniform(.01, 1), n)
	function random_inhibitor(n_activators, n_repressors)
		n_activators == n_repressors ? rand([true, false]) : n_activators < n_repressors
	end
end

struct Gene
	modules::Vector{RegulatoryModule}
	α::Vector{Float64}
	# constructor for JSON3
	Gene(modules::Vector{RegulatoryModule}, α::Vector{<:AbstractFloat}) = new(modules, α)
	function Gene(activators::Vector{<:Integer}, repressors::Vector{<:Integer})
		modules = random_modules(activators, repressors)
		new(modules, random_α(modules))
	end
	
	"""
	activators and repressors: array of indexes for input nodes to distribute to some random number of modules for a single gene.
	The indexes should be for transcription factors only since these are the only ones supposed to affect gene expression directly.
	return: ::Array{RegulatoryModule,1}. Array of regulatory modules for the gene.
	"""
	function random_modules(activators::Vector{<:Integer}, repressors::Vector{<:Integer})::Vector{RegulatoryModule}
		max_n_modules = length(activators) + length(repressors)
		# assign inputs randomly to each module
		module_activators = [[] for _ in 1:max_n_modules]
		module_repressors = [[] for _ in 1:max_n_modules]
		for activator in activators
			push!(module_activators[rand(1:length(module_activators))], activator)
		end
		for repressor in repressors
			push!(module_repressors[rand(1:length(module_repressors))], repressor)
		end

		[RegulatoryModule(module_activators[i], module_repressors[i]) for i in 1:max_n_modules
		if !isempty(module_activators[i]) || !isempty(module_repressors[i])]
	end
	"""
	return: 1D array with length 2^n_modules. activation α values for each unique binding state.
	"""
	function random_α(modules::Vector{RegulatoryModule})::Vector{Float64}
		n_modules = length(modules)
		if n_modules == 0 return [1.] end
		inhibitors = [m.inhibitor for m in modules]
		n_inhibitors = sum(inhibitors)
		# effect on α₀ from each individual module.
		α′ = random_α′(n_modules)
		α′_enhancer = @view α′[.!inhibitors]
		α′_inhibitor = @view α′[inhibitors]
		α′_inhibitor .*= -1  # in-place so the non-view array is modified
		α₀ = random_α₀(n_modules, n_inhibitors)
		
		α_max = α₀ + sum(α′_enhancer)
		α_min = α₀ + sum(α′_inhibitor)
		# make sure that the activation goes at least to 1 in the maximally activated state,
		# when there is at least 1 activator.
		if n_modules > n_inhibitors && α_max < 1
			α′_enhancer[argmin(α′_enhancer)] += 1 - α_max
		end
		# make sure that the activation falls within [0 weak_activation] in the maximally repressed state,
		# if there is a at least one repressor
		if n_inhibitors > 0 && α_min > weak_activation
			# decrease the weakest α′ so that: (α₀ + sum(α′_inhibitor)) in [0 weak_activation]
			α′_inhibitor[argmax(α′_inhibitor)] += random_low_α₀() - α_min
		end
		
		# Set the α for each possible state. State 0 (index 1) will simply become α₀
		α = [α₀ + sum(α′[state]) for state in states(n_modules)]
		clamp!(α, 0, 1)
	end
	
	random_α′(n) = rand(TruncNormal(weak_activation, 1), n)
	"""
	Basal transcription rate.
	"""
	function random_α₀(n_modules, n_inhibitors)
		if n_modules == n_inhibitors return 1  # max expression before inhibitors bind
		elseif n_inhibitors == 0 return random_low_α₀()
		else return random_medium_α₀() end
	end
	random_low_α₀() = rand(TruncatedNormal(0, weak_activation, 0, .05))
	random_medium_α₀() = rand(TruncNormal(weak_activation, 1 - weak_activation))
end

function Base.show(io::IO, m::RegulatoryModule)
	activators = m.inputs[1:m.n_activators]
	repressors = m.inputs[m.n_activators+1:end]
	inputs = []
	if !isempty(activators) push!(inputs, "activators=$activators") end
	if !isempty(repressors) push!(inputs, "repressors=$repressors") end
	print(io, "RegulatoryModule(", join(inputs, ", "), ")")
end
function Base.show(io::IO, g::Gene)
	n_modules = length(g.modules)
	# we take inhibitor status of modules into account when counting activators and repressors of a gene
	activators = [m.inputs[1:m.n_activators] for m in g.modules if !m.inhibitor]
	repressors = [m.inputs[m.n_activators+1:end] for m in g.modules if !m.inhibitor]
	append!(activators, [m.inputs[m.n_activators+1:end] for m in g.modules if m.inhibitor])
	append!(repressors, [m.inputs[1:m.n_activators] for m in g.modules if m.inhibitor])
	activators = vcat(activators...)
	repressors = vcat(repressors...)
	activators = isempty(activators) ? "" : ", activators=$activators"
	repressors = isempty(repressors) ? "" : ", repressors=$repressors"
	print(io, "Gene(n_modules=$n_modules$activators$repressors)")
end


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

"""
Concentration of active protein, which is either phosphorylated or unphosphorylated concentration of the protein, depending on a bool.
"""
ψ(p, ϕ, phos_activation) = @. phos_activation * ϕ + (!phos_activation) * (p-ϕ)

end;

