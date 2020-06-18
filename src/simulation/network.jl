isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")

"""
Network struct containing genes for simulation.
"""

using .ArrayUtils: TruncNormal

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
	ψ₀::Vector{Float64}
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ψ₀::Vector{<:AbstractFloat})
		new(genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ψ₀)
	end
	function Network(genes::Vector{Gene}, Wₚ::Vector{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ψ₀::Vector{<:AbstractFloat})
		Wₚ = reshape(Wₚ, (nₚ+nₜ,nₚ))  # un-flatten matrix
		new(genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ψ₀)
	end
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:Integer})
		λ_phos = random_λ(Wₚ)
		Network(genes, init_Wₚ₊Wₚ₋(genes, Wₚ, λ_phos)..., λ_phos)
	end
	"""
	Make sure to be exact about using either integer or float for Wₚ 
	since a matrix of floats {-1.,0.,1.} will be seen as the exact edge values and not indication of repression, activation, etc.
	"""
	Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}) = Network(genes, Wₚ₊Wₚ₋(Wₚ)..., random_λ(size(Wₚ,1)))
	function Network(genes::Vector{Gene}, Wₚ₊::Matrix{<:AbstractFloat}, Wₚ₋::Matrix{<:AbstractFloat}, λ_phos::Vector)
		nᵥ, nₚ = length(genes), size(Wₚ₊,2)
		nₜ = size(Wₚ₊, 1) - nₚ
		# In the non-dimensionalized model, max_transcription == λ_mRNA and max_translation == λ_prot
		max_transcription = λ_mRNA = random_λ(n)
		max_translation = λ_prot = random_λ(n)
		r₀ = initial_r(max_transcription, λ_mRNA, genes)
		p₀ = initial_p(max_translation, λ_prot, r₀)
		ψ₀ = initial_ψ(Wₚ₊, Wₚ₋, λ_phos, p₀[1:nₜ+nₚ])
		new(genes, Wₚ₊, Wₚ₋, nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ψ₀)
	end
	function Network(Wₜ::Matrix, Wₚ::Matrix)
		nₚ = size(Wₚ,2)
		# == 0 → no edge, > 0 → activator, < 0 → repressor.
		# .+ nₚ because the indexes among all proteins referring to TFs starts after PKs
		genes = [Gene(findall(row .> 0) .+ nₚ, findall(row .< 0) .+ nₚ) for row in eachrow(Wₜ)]
		Network(genes, Wₚ)
	end
	Network(W, nₜ, nₚ) = Network(WₜWₚ(W,nₜ,nₚ)...)
	function Network(net::Network)
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, net.max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ψ₀)
	end
	Base.copy(net::Network) = Network(net)
	"""
	Create a mutant by making a copy of a wildtype network and changing the max transcription level of 1 or more genes.
	"""
	function Network(net::Network, mutate::Integer, value=1e-5)
		max_transcription = copy(net.max_transcription)
		max_transcription[mutate] = value
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
			net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ψ₀)
	end
	function Network(net::Network, mutate::AbstractVector, value=1e-5)
		max_transcription = copy(net.max_transcription)
		mutatable = @view max_transcription[1:net.nₚ+net.nₜ]
		mutatable[mutate] .= value
		new(net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
            net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ψ₀)
	end
	
	random_t½() = TruncNormal(5, 50)
	"""
	We have exponential decay, the half-life and the decay rate are thus related by:
	t½ = ln(2) / λ ⟹
	λ = ln(2) / t½
	"""
	random_λ(n::Int) = log(2) ./ rand(random_t½(), n)
    """
    For λ_phos to have lower values for nodes that are mostly negatively regulated.
    """
    function random_λ(mat::Matrix)
        # weigh by the fraction of regulators that regulate positively.
        positives = sum(mat .> 0; dims=2) |> vec
        negatives = sum(mat .< 0; dims=2) |> vec
        random_λ(size(mat,1)) .* positives ./ (positives .+ negatives)
    end
	
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
	Initial active protein concentrations. Estimated as (p is used in place of ψ)
	0 = dψ/dt = (Wₚ₊ p) (p - ψ) - (Wₚ₋ p + λ_phos) ψ ⟹
	(Wₚ₊ p) p = (Wₚ₊ p - Wₚ₋ p - λ_phos) ψ ⟹
	ψ = Wₚ₊ p p / (Wₚ₊ p - Wₚ₋ p - λ_phos)
	- p: protein concentrations of TFs+PKs
	"""
	function initial_ψ(Wₚ₊::Matrix, Wₚ₋::Matrix, λ_phos::Vector, p::Vector)
		# using ψ = p from the KPs, meaning we use fully active KPs at time=0.
		nₚ = size(Wₚ₊,2)
		Wₚ₊p = Wₚ₊ * p[1:nₚ]
		Wₚ₋p = Wₚ₋ * p[1:nₚ]
		clamp.(Wₚ₊p .* p ./ (Wₚ₊p .- Wₚ₋p .- λ_phos), 0, p)
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
- p: size nᵥ.
- ψₜₚ: size nₜ+nₚ.
"""
drdt(net::Network, r, p, ψₜₚ) = net.max_transcription .* f(net.genes, ψₜₚ) .- net.λ_mRNA .* r
"""
- r: size nᵥ.
- p: size nᵥ.
"""
dpdt(net::Network, r, p) = net.max_translation .* r .- net.λ_prot .* p
"""
- pₜₚ: size nₜ+nₚ.
- ψₜₚ: size nₜ+nₚ.
"""
dψdt(net::Network, pₜₚ, ψₜₚ) = dψdt(net, pₜₚ, ψₜₚ, view(ψₜₚ, 1:net.nₚ))
"""
- pₜₚ: size nₜ+nₚ. Protein concentrations.
- ψₜₚ: size nₜ+nₚ. Active protein concentrations.
- ψₚ: size nₚ. Active protein concentrations.
"""
dψdt(net::Network, pₜₚ, ψₜₚ, ψₚ) = (net.Wₚ₊ * ψₚ) .* (pₜₚ .- ψₜₚ) .- (net.Wₚ₋ * ψₚ .+ net.λ_phos) .* ψₜₚ

