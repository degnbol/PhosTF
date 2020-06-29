isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")

"""
Network struct containing genes for simulation.
"""

using .ArrayUtils: TruncNormal


function init_genes(Wₜ, nₚ)
	# row .== 0 → no edge, row .> 0 → activator, row .< 0 → repressor.
	# .+ nₚ because the indexes among all proteins referring to TFs starts after KPs
	[Gene(findall(row .> 0) .+ nₚ, findall(row .< 0) .+ nₚ) for row in eachrow(Wₜ)]
end


struct Network
	genes::Vector{Gene}
	# KP to TF+KP edges. Using Array instead of Matrix so JSON3 can allow a Vector here before reshape.
	Wₚ₊::Array{Float64}
	Wₚ₋::Array{Float64}
	nᵥ::Int  # number of genes
	nₒ::Int  # number of genes which are not regulators
	nₜ::Int  # number of transcription factors
	nₚ::Int  # number of protein kinases
	max_transcription::Vector{Float64}
	max_translation::Vector{Float64}
	λ_mRNA::Vector{Float64}
	λ_prot::Vector{Float64}
	λ_phos::Vector{Float64} # size nₜ+nₚ
	r₀::Vector{Float64}
	p₀::Vector{Float64}
	a₀::Vector{Float64}
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, a₀::Vector{<:AbstractFloat})
		new(genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, a₀)
	end
	function Network(genes::Vector{Gene}, Wₚ::Vector{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, a₀::Vector{<:AbstractFloat})
		Wₚ = reshape(Wₚ, (nₚ+nₜ,nₚ))  # un-flatten matrix
		new(genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, a₀)
	end
	function Network(genes::Vector{Gene}, Wₚ₊::Matrix{<:AbstractFloat}, Wₚ₋::Matrix{<:AbstractFloat}, λ_phos::Vector)
		nᵥ, nₚ = length(genes), size(Wₚ₊,2)
		nₜ = size(Wₚ₊,1) - nₚ
		# In the non-dimensionalized model, max_transcription == λ_mRNA and max_translation == λ_prot
		max_transcription = λ_mRNA = random_λ(nᵥ)
		max_translation = λ_prot = random_λ(nᵥ)
		r₀ = initial_r(max_transcription, λ_mRNA, genes)
		p₀ = initial_p(max_translation, λ_prot, r₀)
		a₀ = initial_ψ(Wₚ₊, Wₚ₋, λ_phos, p₀[1:nₜ+nₚ])
		new(genes, Wₚ₊, Wₚ₋, nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, a₀)
	end
	"""
	Make sure to be exact about using either integer or float for Wₚ 
	since a matrix of floats {-1.,0.,1.} will be seen as the exact edge values and not indication of repression, activation, etc.
	"""
	Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}) = Network(genes, Wₚ₊Wₚ₋(Wₚ)..., random_λ(size(Wₚ,1)))
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:Integer})
		λ_phos = random_λ(Wₚ)
		Network(genes, init_Wₚ₊Wₚ₋(genes, Wₚ, λ_phos)..., λ_phos)
	end
	Network(Wₜ::Matrix, Wₚ::Matrix) = Network(init_genes(Wₜ, size(Wₚ,2)), Wₚ)
	Network(W, nₜ::Integer, nₚ::Integer) = Network(WₜWₚ(W,nₜ,nₚ)...)
	function Network(net::Network)
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, net.max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.a₀)
	end
	"""
	Create a mutant by making a copy of a wildtype network and changing the max transcription level of 1 or more genes.
	"""
	function Network(net::Network, mutate::Integer, value=1e-5)
		max_transcription = copy(net.max_transcription)
		max_transcription[mutate] = value
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
			net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.a₀)
	end
	function Network(net::Network, mutate::AbstractVector, value=1e-5)
		max_transcription = copy(net.max_transcription)
		mutatable = @view max_transcription[1:net.nₚ+net.nₜ]
		mutatable[mutate] .= value
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
            net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.a₀)
	end
	Base.copy(net::Network) = Network(net)
	
	"""
	Initial mRNA. Estimated as
	0 = dr/dt = max_transcription*f(ψ) -λ_mRNA*r ⟹ r = m*f(ψ) / λ_mRNA  
	where we use p = 0 ⟹ ψ = 0  
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
	Initial protein activities.
	- p: protein concentrations
	"""
	function initial_a(Wₚ₊::Matrix, Wₚ₋::Matrix, λ_phos::Vector, p::Vector)
		nₚ = size(Wₚ₊,2)
        nₚ₊= sum(Wₚ₊ .> 0; dims=2) |> vec
        nₚ₋= sum(Wₚ₋ .> 0; dims=2) |> vec
		# Sum the effects from regulators when using ψ = p from the KPs, which is fully active KPs.
		Wₚ₊p = Wₚ₊ * p[1:nₚ]
		Wₚ₋p = Wₚ₋ * p[1:nₚ]
        # Fraction of proteins that are activated on a scale [0,1] as weighted average.
        # We subtract deactivations from max activity level, and add activations to min activity level.
        a = @. (Wₚ₊p * nₚ₊ + (1 - Wₚ₋p) * nₚ₋) / (nₚ₊ + nₚ₋)
        # we get NaN if a protein is not regulated. In that case simply let it be fully active.
        a[isnan.(a)] .= 1
        a
	end
end

Base.show(io::IO, net::Network) = print(io, "Network(nᵥ=$(net.nᵥ), nₜ=$(net.nₜ), nₚ=$(net.nₚ))")

function Wₚ₊Wₚ₋(Wₚ::AbstractMatrix)
	Wₚ₊ = +Wₚ; Wₚ₊[Wₚ₊ .<= 0] .= 0
	Wₚ₋ = -Wₚ; Wₚ₋[Wₚ₋ .<= 0] .= 0
	Wₚ₊, Wₚ₋
end

"""
- r: size nᵥ.
- pₚₜ: size nₚ+nₜ.
- aₚₜ: size nₚ+nₜ.
"""
drdt(net::Network, r, pₚₜ, aₚₜ) = net.max_transcription .* f(net.genes, aₚₜ .* pₚₜ) .- net.λ_mRNA .* r
"""
- r: size nᵥ.
- p: size nᵥ.
"""
dpdt(net::Network, r, p) = net.max_translation .* r .- net.λ_prot .* p
"""
- pₚₜ: size nₚ+nₜ.
- ψₚₜ: size nₚ+nₜ.
"""
dadt(net::Network, pₚₜ, ψₚₜ) = dψdt(net, pₜₚ, ψₜₚ, view(ψₜₚ, 1:net.nₚ)) ./ pₚₜ
"""
- pₜₚ: size nₜ+nₚ. Protein concentrations.
- ψₜₚ: size nₜ+nₚ. Active protein concentrations.
- ψₚ: size nₚ. Active protein concentrations.
"""
dψdt(net::Network, pₜₚ, ψₜₚ, ψₚ) = (net.Wₚ₊ * ψₚ) .* (pₜₚ .- ψₜₚ) .- (net.Wₚ₋ * ψₚ .+ net.λ_phos) .* ψₜₚ

