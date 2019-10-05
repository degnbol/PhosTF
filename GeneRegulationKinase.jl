#!/usr/bin/env julia
"""
Backup version which assumes all phosphorylation regulators are kinases, so no negative (phosphatase) edges.
Module intended to hold all structs with data defining our gene regulation network and its regulation mechanisms. All functions for initializing the values are found here, including randomized initialization.
"""
module GeneRegulation
using Distributions: Normal, Uniform, mean, TruncatedNormal
import JSON3
include("utilities/ArrayUtils.jl"); using .ArrayUtils
include("Model.jl"); using .Model: WₜWₚ

export Network
export drdt, dpdt, dϕdt
export mutate

const weak_activation = .25
const noise_activation = .25weak_activation

struct RegulatoryModule
	n_activators::Int
	n_repressors::Int
	inputs::Vector{Int}  # int indexes for activators then repressors among all proteins
	ν::Vector{Float64}
	k::Vector{Float64}
	complex::Bool
	inhibitor::Bool
	# constructor for JSON3
	function RegulatoryModule(n_activators::Int, n_repressors::Int, inputs::Vector{Int}, ν::Vector{Float64}, k::Vector{Float64}, complex::Bool, inhibitor::Bool)
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
	random_ν(n::Int) = rand(TruncatedNormal(2, 2, 1, 10), n)
	"""
	Dissociation constant k for each input protein.
	"""
	random_k(n::Int) = rand(Uniform(.01, 1), n)
	function random_inhibitor(n_activators, n_repressors)
		n_activators == n_repressors ? rand([true, false]) : n_activators < n_repressors
	end
	
	function Base.show(io::IO, m::RegulatoryModule)
		activators = m.inputs[1:m.n_activators]
		repressors = m.inputs[m.n_activators+1:end]
		inputs = []
		if !isempty(activators) push!(inputs, "activators=$activators") end
		if !isempty(repressors) push!(inputs, "repressors=$repressors") end
		print("RegulatoryModule(", join(inputs, ", "), ")")
	end
end

struct Gene
	modules::Vector{RegulatoryModule}
	α::Vector{Float64}
	# constructor for JSON3
	Gene(modules::Vector{RegulatoryModule}, α::Vector{Float64}) = new(modules, α)
	function Gene(activators::Vector{Int}, repressors::Vector{Int})
		modules = random_modules(activators, repressors)
		new(modules, random_α(modules))
	end
	
	"""
	activators and repressors: array of indexes for input nodes to distribute to some random number of modules for a single gene.
	The indexes should be for transcription factors only since these are the only ones supposed to affect gene expression directly.
	return: ::Array{RegulatoryModule,1}. Array of regulatory modules for the gene.
	"""
	function random_modules(activators::Vector{Int}, repressors::Vector{Int})::Vector{RegulatoryModule}
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
		α′_inhibitor[:] *= -1  # in-place so the main array is modified
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
		print("Gene(n_modules=$n_modules$activators$repressors)")
	end
end

struct Network
	genes::Vector{Gene}
	# PK to TF+PK edges. Using Array instead of Matrix so JSON3 can allow a Vector here before reshape.
	Wₚ::Array{Float64}
	n::Int  # number of genes
	nₜ::Int  # number of transcription factors
	nₚ::Int  # number of protein kinases
	max_transcription::Vector{Float64}
	max_translation::Vector{Float64}
	λ_mRNA::Vector{Float64}
	λ_prot::Vector{Float64}
	λ_phos::Vector{Float64}
	r₀::Vector{Float64}
	p₀::Vector{Float64}
	ϕ₀::Vector{Float64}
	phos_activation::BitVector
	function Network(genes::Vector{Gene}, Wₚ::Matrix{Float64}, n::Int, nₜ::Int, nₚ::Int, max_transcription::Vector{Float64}, max_translation::Vector{Float64}, λ_mRNA::Vector{Float64}, λ_prot::Vector{Float64}, λ_phos::Vector{Float64}, r₀::Vector{Float64}, p₀::Vector{Float64}, ϕ₀::Vector{Float64}, phos_activation::BitVector)
		new(genes, Wₚ, n, nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
	end
	function Network(genes::Vector{Gene}, Wₚ::Vector{Float64}, n::Int, nₜ::Int, nₚ::Int, max_transcription::Vector{Float64}, max_translation::Vector{Float64}, λ_mRNA::Vector{Float64}, λ_prot::Vector{Float64}, λ_phos::Vector{Float64}, r₀::Vector{Float64}, p₀::Vector{Float64}, ϕ₀::Vector{Float64}, phos_activation::BitVector)
		Wₚ = reshape(Wₚ, (nₚ+nₜ,nₚ))  # un-flatten matrix
		new(genes, Wₚ, n, nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
	end
	function Network(genes::Vector{Gene}, Wₚ::Matrix{Float64})
		n, nₚ = length(genes), size(Wₚ, 2)
		nₜ = size(Wₚ, 1) - nₚ
		@assert n >= nₜ + nₚ
		# In the non-dimensionalized model, max_transcription == λ_mRNA and max_translation == λ_prot
		max_transcription = λ_mRNA = random_λ(n)
		max_translation = λ_prot = random_λ(n)
		λ_phos = random_λ(nₚ+nₜ)
		r₀ = initial_r(max_transcription, λ_mRNA, genes)
		p₀ = initial_p(max_translation, λ_prot, r₀)
		ϕ₀ = initial_ϕ(Wₚ, λ_phos, p₀[1:nₚ+nₜ])
		# activated by phosphorylation if there is kinase regulation on a protein (and more kinases than phosphatases)
		phos_activation = [vec(sum(Wₚ, dims=2)) .> 0; falses(n-(nₚ+nₜ))]
		new(genes, Wₚ, n, nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
	end
	function Network(Wₜ::Matrix, Wₚ::Matrix{Float64})
		nₚ = size(Wₚ,2)
		# == 0 → no edge, > 0 → activator, < 0 → repressor.
		# .+ nₚ because the indexes among all proteins referring to TFs starts after PKs
		genes = [Gene(findall(row .> 0) .+ nₚ, findall(row .< 0) .+ nₚ) for row in eachrow(Wₜ)]
		Network(genes, Wₚ)
	end
	Network(W, nₜ, nₚ) = Network(WₜWₚ(W,nₜ,nₚ)...)
	function Network(net::Network)
		new(net.genes, net.Wₚ, net.n, net.nₜ, net.nₚ, net.max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	Base.copy(net::Network) = Network(net)
	"""
	Create a mutant by making a copy of a wildtype network and changing the max transcription level of 1 or more genes.
	"""
	function Network(net::Network, mutate::Int, value=1e-7)
		max_transcription = copy(net.max_transcription)
		max_transcription[mutate] = value
		new(net.genes, net.Wₚ, net.n, net.nₜ, net.nₚ, max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	function Network(net::Network, mutate::T, value=1e-7) where T<:AbstractVector
		max_transcription = copy(net.max_transcription)
		mutatable = @view max_transcription[1:net.nₚ+net.nₜ]
		mutatable[mutate] .= value
		new(net.genes, net.Wₚ, net.n, net.nₜ, net.nₚ, max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	
	"""
	We have exponential decay, the half-life and the decay rate are thus related by:
	t½ = ln(2) / λ ⟹
	λ = ln(2) / t½
	"""
	random_λ(n) = log(2) ./ rand(random_t½(), n)
	random_t½() = TruncNormal(5, 50)
	
	"""
	Initial mRNA. Estimated as
	0 = dr/dt = max_transcription*f(ϕ) -λ_mRNA*r ⟹ r = m*f(ϕ) / λ_mRNA
	where we use p = 0 ⟹ ϕ = 0
	"""
	function initial_r(max_transcription::Vector, λ_mRNA::Vector, genes::Vector{Gene})
		max_transcription .* f(genes, zeros(length(genes))) ./ λ_mRNA
	end
	"""
	Initial protein concentrations. Estimated as
	0 = dp/dt = max_translation*r - λ_prot*p ⟹ p = max_translation*r / λ_prot
	"""
	initial_p(max_translation::Vector, λ_prot::Vector, r::Vector) = max_translation .* r ./ λ_prot
	"""
	Initial phosphorylated protein concentrations. Estimated as (p is used in place of ϕ)
	0 = dϕ/dt = (Wₚ * p) * (p - ϕ) - λ_phos * ϕ ⟹
	(Wₚ * p - λ_phos) * ϕ = (Wₚ * p) * p
	ϕ = (Wₚ * 1) * p / (Wₚ * 1 - λ_phos)
	p: protein concentrations of TFs+PKs
	"""
	function initial_ϕ(Wₚ::Matrix, λ_phos::Vector, p::Vector)
		# using ϕ = p from the PKs, meaning we use fully phosphorylated PKs as the origin of calculation
		nₚ = size(Wₚ,2)
		Wp = Wₚ * (p[1:nₚ])
		clamp.(Wp .* p ./ (Wp .- λ_phos), 0, p)
	end
	
	Base.show(io::IO, net::Network) = print("Network(n=$(net.n), nₜ=$(net.nₜ), nₚ=$(net.nₚ))")
end


"""
State s is interpreted as a binary number, where bit k indicates whether module k is active (True) or inactive (False) in this state.
n: number of modules to combine in each unique way
"""
states(n) = (binary(i-1, n) for i in 1:2^n)

"""
Mean activation μ given ψ.
ψ: 1D array. Active nondim concentrations.
"""
function μ(m::RegulatoryModule, ψ::Vector{Float64})
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
μ(ms::Vector{RegulatoryModule}, ψ::Vector{Float64}) = [μ(m, ψ) for m in ms]
"""
Fraction of max activation for a given gene when active TFs are found at a given concentration.
"""
function f(gene::Gene, ψ::Vector{Float64})
	if isempty(gene.modules) return 1 end
	μs = μ(gene.modules, ψ)
	# get P{state} for all states, each state is a unique combination of modules
	P = [prod(μs[state]) * prod(1 .- μs[.!state]) for state in states(length(gene.modules))]
	sum(gene.α .* P)
end
f(genes::Vector{Gene}, ψ::Vector{Float64}) = [f(gene, ψ) for gene in genes]

"""
Concentration of active protein, which is either phosphorylated or unphosphorylated concentration of the protein, depending on a bool.
"""
ψ(p, ϕ, phos_activation) = @. phos_activation * ϕ + (!phos_activation) * (p-ϕ)

drdt(network::Network, r, p, ϕ) = network.max_transcription .* f(network.genes, ψ(p, ϕ, network.phos_activation)) .- network.λ_mRNA .* r
dpdt(network::Network, r, p) = network.max_translation .* r .- network.λ_prot .* p
dϕdt(network::Network, p, ϕ) = (network.Wₚ * ϕ[1:network.nₚ]) .* (p .- ϕ) .- network.λ_phos .* ϕ

"Estimate the effect on f for all genes when a given TF has either ϕ=weak or ϕ=strong."
function estimate_Wₜ(net::Network, i::Int, basal_activation::Float64)
	basal = fill(basal_activation, net.n)
	ψ = copy(basal); ψ[i] = 1
	f(net.genes, ψ) - f(net.genes, basal)
end
function estimate_Wₜ(net::Network, i::Int)
	mean(estimate_Wₜ(net, i, activation) for activation in [0, noise_activation, weak_activation])
end
"Estimate Wₜ by comparing the effect on f when any TF has ϕ=weak or ϕ=strong."
estimate_Wₜ(net::Network) = hcat([estimate_Wₜ(net, i) for i in net.nₚ+1:net.nₚ+net.nₜ]...)

# define struct types for JSON3 to be able to read/write them
JSON3.StructType(::Type{RegulatoryModule}) = JSON3.Struct()
JSON3.StructType(::Type{Gene}) = JSON3.Struct()
JSON3.StructType(::Type{Network}) = JSON3.Struct()

end;

