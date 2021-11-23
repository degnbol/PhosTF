#!/usr/bin/env julia
include("../utilities/ArrayUtils.jl") # different scope, have to re-include
using Distributions: TruncatedNormal
using .ArrayUtils: TruncNormal, binary

"""
State s is interpreted as a binary number, where bit k indicates whether module k is active (true) or inactive (false) in this state. The first returned value is all falses and the last is all trues.
n: number of modules to combine in each unique way
"""
states(n) = (binary(i-1, n) for i in 1:2^n)

"""
Gene struct for gene regulation simulation.
"""
struct Gene
	modules::Vector{RegulatoryModule} # regulatory modules regulating this gene
	α::Vector{Float64} # values associated with each state of regulation
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
		# α′ is the effect added on top of α₀ by each individual module.
		α′ = random_α′(n_modules)
		α′_enhancer = @view α′[.!inhibitors]
		α′_inhibitor = @view α′[inhibitors]
		α′_inhibitor .*= -1  # in-place so the non-view array is modified
		α₀ = random_α₀(n_modules, n_inhibitors)
		
		α_max = α₀ + sum(α′_enhancer)
		α_min = α₀ + sum(α′_inhibitor)
		# if there are ≥1 activators then
		# make sure that activation α ≥ random_high_α in the maximally activated state. In GNW it was 1 instead of random high α
		if n_modules > n_inhibitors && α_max < 1
			α′_enhancer[argmin(α′_enhancer)] += random_high_α₀() - α_max
		end
		# if there are ≥1 repressors then
		# make sure that activation α ∈ [0,weak_activation] in the maximally repressed state.
		if n_inhibitors > 0 && α_min > weak_activation
			# decrease the weakest α′ so that: (α₀ + sum(α′_inhibitor)) ∈ [0,weak_activation]
			α′_inhibitor[argmax(α′_inhibitor)] += random_low_α₀() - α_min
		end
		
		# Set the α for each possible state. State 0 (index 1) will simply become α₀
		α = [α₀ + sum(α′[state]) for state in states(n_modules)]
		clamp!(α, 0, 1)
	end
	
	random_α′(n) = rand(TruncNormal(weak_activation, 1.), n)
	"""
	Basal transcription rate.
	"""
	function random_α₀(n_modules, n_inhibitors)
		if n_modules == n_inhibitors return 1.  # max expression before inhibitors bind
		elseif n_inhibitors == 0 return random_low_α₀()
		else return random_medium_α₀() end
	end
	"GNW has both used (0,.05,0,weak) and (.001,.05,.001,weak)"
    random_low_α₀()    = rand(truncated(Normal(.05, .05), .01, weak_activation))
	random_medium_α₀() = rand(TruncNormal(weak_activation, 1. - weak_activation))
	random_high_α₀()   = rand(TruncNormal(strong_activation, 1.))
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
