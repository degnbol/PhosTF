isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
using Distributions: Uniform, truncated, Normal

"""
Regulatory module struct for gene regulation simulation.
"""
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
	GNW defaults are μ=2, σ=2, min=1, max=10
	"""
    random_ν(n::Integer) = rand(truncated(Normal(2., 2.), .5, 8.), n)
	"""
	Dissociation constant k for each input protein.
	"""
	random_k(n::Integer) = rand(Uniform(.05, .7), n) # GNW uses (.01, 1.)
	function random_inhibitor(n_activators, n_repressors)
		n_activators == n_repressors ? rand([true, false]) : n_activators < n_repressors
	end
end

function Base.show(io::IO, m::RegulatoryModule)
	activators = m.inputs[1:m.n_activators]
	repressors = m.inputs[m.n_activators+1:end]
	inputs = []
	if !isempty(activators) push!(inputs, "activators=$activators") end
	if !isempty(repressors) push!(inputs, "repressors=$repressors") end
	print(io, "RegulatoryModule(", join(inputs, ", "), ")")
end
