if !isdefined(Main, :ArrayUtils) include("../utilities/ArrayUtils.jl") end

"""
Network struct containing genes for simulation.
"""

using .ArrayUtils: TruncNormal

struct Network
	genes::Vector{Gene}
	# PK to TF+PK edges. Using Array instead of Matrix so JSON3 can allow a Vector here before reshape.
	Wₚₖ::Array{Float64}
	Wₚₚ::Array{Float64}
	n::Int  # number of genes
	nₓ::Int  # number of genes which are not regulators
	nₜ::Int  # number of transcription factors
	nₚ::Int  # number of protein kinases
	max_transcription::Vector{Float64}
	max_translation::Vector{Float64}
	λ_mRNA::Vector{Float64}
	λ_prot::Vector{Float64}
	λ_phos::Vector{Float64} # size nₜ+nₚ
	r₀::Vector{Float64}
	p₀::Vector{Float64}
	ϕ₀::Vector{Float64}
	phos_activation::BitVector # size nₜ+nₚ
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}, n::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ϕ₀::Vector{<:AbstractFloat}, phos_activation::BitVector)
		new(genes, WₚₖWₚₚ(Wₚ)..., n, n-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
	end
	function Network(genes::Vector{Gene}, Wₚ::Vector{<:AbstractFloat}, n::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ_phos::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ϕ₀::Vector{<:AbstractFloat}, phos_activation::BitVector)
		Wₚ = reshape(Wₚ, (nₚ+nₜ,nₚ))  # un-flatten matrix
		new(genes, WₚₖWₚₚ(Wₚ)..., n, n-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
	end
	function Network(genes::Vector{Gene}, Wₚ::Matrix{<:Integer})
		phos_activation=vec(init_phos_activation(Wₚ))
		λ_phos = random_λ(size(Wₚ,1))
		Network(genes, init_phos_edges(genes, Wₚ, phos_activation, λ_phos)..., phos_activation, λ_phos)
	end
	"""
	Make sure to be exact about using either integer or float for Wₚ 
	since a matrix of floats {-1.,0.,1.} will be seen as the exact edge values and not indication of repression, activation, etc.
	"""
	Network(genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}) = Network(genes, WₚₖWₚₚ(Wₚ)..., vec(init_phos_activation(Wₚ)), random_λ(size(Wₚ,1)))
	function Network(genes::Vector{Gene}, Wₚₖ::Matrix{<:AbstractFloat}, Wₚₚ::Matrix{<:AbstractFloat}, phos_activation::BitVector, λ_phos::Vector)
		n, nₚ = length(genes), size(Wₚₖ,2)
		nₜ = size(Wₚₖ, 1) - nₚ
		# In the non-dimensionalized model, max_transcription == λ_mRNA and max_translation == λ_prot
		max_transcription = λ_mRNA = random_λ(n)
		max_translation = λ_prot = random_λ(n)
		r₀ = initial_r(max_transcription, λ_mRNA, genes)
		p₀ = initial_p(max_translation, λ_prot, r₀)
		ϕ₀ = initial_ϕ(Wₚₖ, Wₚₚ, λ_phos, p₀[1:nₜ+nₚ])
		new(genes, Wₚₖ, Wₚₚ, n, n-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ_phos, r₀, p₀, ϕ₀, phos_activation)
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
		new(net.genes, net.Wₚₖ, net.Wₚₚ, net.n, net.n-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, net.max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	Base.copy(net::Network) = Network(net)
	"""
	Create a mutant by making a copy of a wildtype network and changing the max transcription level of 1 or more genes.
	"""
	function Network(net::Network, mutate::Integer, value=1e-5)
		max_transcription = copy(net.max_transcription)
		max_transcription[mutate] = value
		new(net.genes, net.Wₚₖ, net.Wₚₚ, net.n, net.n-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
			net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	function Network(net::Network, mutate::AbstractVector, value=1e-5)
		max_transcription = copy(net.max_transcription)
		mutatable = @view max_transcription[1:net.nₚ+net.nₜ]
		mutatable[mutate] .= value
		new(net.genes, net.Wₚₖ, net.Wₚₚ, net.n, net.n-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
			net.λ_mRNA, net.λ_prot, net.λ_phos, net.r₀, net.p₀, net.ϕ₀, net.phos_activation)
	end
	
	"""
	We have exponential decay, the half-life and the decay rate are thus related by:
	t½ = ln(2) / λ ⟹
	λ = ln(2) / t½
	"""
	random_λ(n) = log(2) ./ rand(random_t½(), n)
	random_t½() = TruncNormal(5, 50)
	
	"Activated by phosphorylation if there is phos regulation on a protein (and at least as many kinases as phosphatases)"
	init_phos_activation(Wₚ::Matrix) = (sum(abs.(Wₚ), dims=2) .> 0) .& (sum(Wₚ, dims=2) .>= 0)
	
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
	0 = dϕ/dt = (Wₚₖ p) (p - ϕ) - (Wₚₚ p + λ_phos) ϕ ⟹
	(Wₚₖ p) p = (Wₚₖ p - Wₚₚ p - λ_phos) ϕ ⟹
	ϕ = Wₚₖ p p / (Wₚₖ p - Wₚₚ p - λ_phos)
	- p: protein concentrations of TFs+PKs
	"""
	function initial_ϕ(Wₚₖ::Matrix, Wₚₚ::Matrix, λ_phos::Vector, p::Vector)
		# using ϕ = p from the PKs, meaning we use fully phosphorylated PKs as the origin of calculation
		nₚ = size(Wₚₖ,2)
		Wₚₖp = Wₚₖ * (p[1:nₚ])
		Wₚₚp = Wₚₚ * (p[1:nₚ])
		clamp.(Wₚₖp .* p ./ (Wₚₖp .- Wₚₚp .- λ_phos), 0, p)
	end
end

Base.show(io::IO, net::Network) = print(io, "Network(n=$(net.n), nₜ=$(net.nₜ), nₚ=$(net.nₚ))")

function WₚₖWₚₚ(Wₚ::AbstractMatrix)
	Wₚₖ = +Wₚ; Wₚₖ[Wₚₖ .<= 0] .= 0
	Wₚₚ = -Wₚ; Wₚₚ[Wₚₚ .<= 0] .= 0
	Wₚₖ, Wₚₚ
end

"""
- r: size n.
- p: size n.
- ϕₜₚ: size nₜ+nₚ.
"""
function drdt(net::Network, r, p, ϕₜₚ)
	net.max_transcription .* f(net.genes, ψ(view(p,1:net.nₜ+net.nₚ), ϕₜₚ, net.phos_activation)) .- net.λ_mRNA .* r
end
"""
- r: size n.
- p: size n.
"""
dpdt(net::Network, r, p) = net.max_translation .* r .- net.λ_prot .* p
"""
- pₜₚ: size nₜ+nₚ.
- ϕₜₚ: size nₜ+nₚ.
"""
dϕdt(net::Network, pₜₚ, ϕₜₚ) = dϕdt(net, pₜₚ, ϕₜₚ, view(ψ(pₜₚ, ϕₜₚ, net.phos_activation), 1:net.nₚ))
"""
- pₜₚ: size nₜ+nₚ. Protein concentrations.
- ϕₜₚ: size nₜ+nₚ. Phosphorylated protein concentrations.
- ψₚ: size nₚ. Active protein concentrations.
"""
dϕdt(net::Network, pₜₚ, ϕₜₚ, ψₚ) = (net.Wₚₖ * ψₚ) .* (pₜₚ .- ϕₜₚ) .- (net.Wₚₚ * ψₚ .+ net.λ_phos) .* ϕₜₚ

