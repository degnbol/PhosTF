#!/usr/bin/env julia
include("Gene.jl")
Main.@use "utilities/ArrayUtils"
Main.@use "inference/Model"
using Main.Model: nₜnₚnₒ

"""
Network struct containing genes for simulation.
"""

random_t½() = Main.ArrayUtils.TruncNormal(5, 50)
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

include("phos_edges.jl") # provides init_Wₚ₊Wₚ₋ and uses random_λ

function init_genes(Wₜ)
    # Each row in Wₜ is edges from each T to a single node. The nodes in rows should be in order T, P, O.
	# row .== 0 → no edge, row .> 0 → activator, row .< 0 → repressor.
	# .+ nₚ because the indexes among all proteins referring to Ts starts after Ps
	[Gene(findall(row .> 0), findall(row .< 0)) for row in eachrow(Wₜ)]
end

default_names(nₜ::Integer, nₚ::Integer, nₒ::Integer) = [["T$i" for i in 1:nₜ]; ["P$i" for i in 1:nₚ]; ["O$i" for i in 1:nₒ]]

function Wₚ₊Wₚ₋(Wₚ::AbstractMatrix)
	Wₚ₊ = +Wₚ; Wₚ₊[Wₚ₊ .<= 0] .= 0
	Wₚ₋ = -Wₚ; Wₚ₋[Wₚ₋ .<= 0] .= 0
	Wₚ₊, Wₚ₋
end


struct Network
    # name of each gene/node
    names::Vector{String}
    # containing information about how each gene (every node) is bound and regulated by Fs.
	genes::Vector{Gene}
	# P to R edges. Using Array instead of Matrix so JSON3 can allow a Vector here before reshape.
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
	λ₊::Vector{Float64} # size nₜ+nₚ
	λ₋::Vector{Float64} # size nₜ+nₚ
	r₀::Vector{Float64}
	p₀::Vector{Float64}
	ψ₀::Vector{Float64}
    function Network(names::Vector{String}, genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ₊::Vector{<:AbstractFloat}, λ₋::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ψ₀::Vector{<:AbstractFloat})
		new(names, genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ₊, λ₋, r₀, p₀, ψ₀)
	end
	function Network(names::Vector{String}, genes::Vector{Gene}, Wₚ::Vector{<:AbstractFloat}, nᵥ::Integer, nₜ::Integer, nₚ::Integer, 
		max_transcription::Vector{<:AbstractFloat}, max_translation::Vector{<:AbstractFloat}, 
		λ_mRNA::Vector{<:AbstractFloat}, λ_prot::Vector{<:AbstractFloat}, λ₊::Vector{<:AbstractFloat}, λ₋::Vector{<:AbstractFloat}, 
		r₀::Vector{<:AbstractFloat}, p₀::Vector{<:AbstractFloat}, ψ₀::Vector{<:AbstractFloat})
		Wₚ = reshape(Wₚ, (nₜ+nₚ, nₚ))  # un-flatten matrix
		new(names, genes, Wₚ₊Wₚ₋(Wₚ)..., nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ₊, λ₋, r₀, p₀, ψ₀)
	end
	function Network(names::Vector{String}, genes::Vector{Gene}, Wₚ₊::Matrix{<:AbstractFloat}, Wₚ₋::Matrix{<:AbstractFloat}, λ₊::Vector, λ₋::Vector)
		nᵥ, nₚ = length(genes), size(Wₚ₊, 2)
		nₜ = size(Wₚ₊, 1) - nₚ
		# In the non-dimensionalized model, max_transcription == λ_mRNA and max_translation == λ_prot
		max_transcription = λ_mRNA = random_λ(nᵥ)
		max_translation = λ_prot = random_λ(nᵥ)
		r₀ = initial_r(max_transcription, λ_mRNA, genes)
		p₀ = initial_p(max_translation, λ_prot, r₀)
		ψ₀ = initial_ψ(Wₚ₊, Wₚ₋, p₀[1:nₜ+nₚ])
		new(names, genes, Wₚ₊, Wₚ₋, nᵥ, nᵥ-(nₜ+nₚ), nₜ, nₚ, max_transcription, max_translation, λ_mRNA, λ_prot, λ₊, λ₋, r₀, p₀, ψ₀)
	end
	"""
	Make sure to be exact about using either integer or float for Wₚ 
	since a matrix of floats {-1.,0.,1.} will be seen as the exact edge values and not indication of repression, activation, etc.
	"""
	Network(names::Vector{String}, genes::Vector{Gene}, Wₚ::Matrix{<:AbstractFloat}) = Network(names, genes, Wₚ₊Wₚ₋(Wₚ)..., random_λ(Wₚ)...)
	function Network(names, genes::Vector{Gene}, Wₚ::Matrix{<:Integer})
		λ₊, λ₋ = random_λ(Wₚ)
		Network(names, genes, init_Wₚ₊Wₚ₋(genes, Wₚ, λ₊, λ₋)..., λ₊, λ₋)
	end
    Network(Wₜ::Matrix, Wₚ::Matrix; names=default_names(nₜnₚnₒ(Wₜ, Wₚ)...)) = Network(names, init_genes(Wₜ), Wₚ)
    Network(W, nₜ::Integer, nₚ::Integer; names=default_names(nₜ, nₚ, size(W,1)-(nₜ+nₚ))) = Network(WₜWₚ(W, nₜ, nₚ)...; names=names)
	function Network(net::Network)
		new(net.names, net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, net.max_transcription, net.max_translation, net.λ_mRNA, net.λ_prot, net.λ₊, net.λ₋, net.r₀, net.p₀, net.ψ₀)
	end
	"""
	Create a mutant by making a copy of a wildtype network and changing the max transcription level of 1 or more genes.
	"""
	function Network(net::Network, mutate::Integer, value=1e-5)
		max_transcription = copy(net.max_transcription)
		max_transcription[mutate] = value
		new(net.names, net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
			net.λ_mRNA, net.λ_prot, net.λ₊, net.λ₋, net.r₀, net.p₀, net.ψ₀)
	end
	function Network(net::Network, mutate::AbstractVector, value=1e-5)
		max_transcription = copy(net.max_transcription)
		mutatable = @view max_transcription[1:net.nₜ+net.nₚ]
		mutatable[mutate] .= value
		new(net.names, net.genes, net.Wₚ₊, net.Wₚ₋, net.nᵥ, net.nᵥ-(net.nₜ+net.nₚ), net.nₜ, net.nₚ, max_transcription, net.max_translation, 
            net.λ_mRNA, net.λ_prot, net.λ₊, net.λ₋, net.r₀, net.p₀, net.ψ₀)
	end
	Base.copy(net::Network) = Network(net)
	
	"""
	Initial mRNA. Estimated as
	0 = dr/dt = max_transcription*f(ψ) -λ_mRNA*r ⟹  r = m*f(ψ) / λ_mRNA  
	where we use p = 0 ⟹ ψ = 0  
	"""
	function initial_r(max_transcription::Vector, λ_mRNA::Vector, genes::Vector{Gene})
        # f(gene, ψ::Vector) is a function from Gene for its fraction of max activation.
        max_transcription .* f.(genes, Ref(zeros(length(genes)))) ./ λ_mRNA
	end
	"""
	Initial protein concentrations. Estimated as
	0 = dp/dt = max_translation*r - λ_prot*p ⟹  p = max_translation*r / λ_prot  
	"""
	initial_p(max_translation::Vector, λ_prot::Vector, r::Vector) = max_translation .* r ./ λ_prot
	"""
	Initial active protein concentrations.
	- p: protein concentrations of R
	"""
	function initial_ψ(Wₚ₊::Matrix, Wₚ₋::Matrix, p::Vector)
		nₚ = size(Wₚ₊, 2)
        nₜ = size(Wₚ₊, 1) - nₚ
        nₚ₊= sum(Wₚ₊ .> 0; dims=2) |> vec
        nₚ₋= sum(Wₚ₋ .> 0; dims=2) |> vec
		# Sum the effects from regulators when using ψ = p from the Ps, which is fully active Ps.
		Wₚ₊p = Wₚ₊ * p[nₜ+1:nₜ+nₚ]
		Wₚ₋p = Wₚ₋ * p[nₜ+1:nₜ+nₚ]
        # Fraction of proteins that are activated on a scale [0,1] as weighted average.
        # We subtract deactivations from max activity level, and add activations to min activity level.
        a = @. (Wₚ₊p * nₚ₊ + (1 - Wₚ₋p) * nₚ₋) / (nₚ₊ + nₚ₋)
        # we get NaN if a protein is not regulated. In that case simply let it be fully active.
        a[isnan.(a)] .= 1
        a .* p 
	end
end

Base.show(io::IO, net::Network) = print(io, "Network(nᵥ=$(net.nᵥ), nₜ=$(net.nₜ), nₚ=$(net.nₚ))")


"""
- r: size nᵥ.
- p: size nᵥ.
- ψₜₚ: size nₜ+nₚ.
"""
drdt(net::Network, r, p, ψₜₚ) = net.max_transcription .* f.(net.genes, Ref(ψₜₚ)) .- net.λ_mRNA .* r
"""
- r: size nᵥ.
- p: size nᵥ.
"""
dpdt(net::Network, r, p) = net.max_translation .* r .- net.λ_prot .* p
"""
- pₜₚ: size nₜ+nₚ.
- ψₜₚ: size nₜ+nₚ.
"""
dψdt(net::Network, pₜₚ, ψₜₚ) = dψdt(net, pₜₚ, ψₜₚ, view(ψₜₚ, net.nₜ+1:net.nₜ+net.nₚ))
"""
- pₜₚ: size nₜ+nₚ. Protein concentrations.
- ψₜₚ: size nₜ+nₚ. Active protein concentrations.
- ψₚ: size nₚ. Active protein concentrations.
"""
dψdt(net::Network, pₜₚ, ψₜₚ, ψₚ) = (net.Wₚ₊ * ψₚ .+ net.λ₊) .* (pₜₚ .- ψₜₚ) .- (net.Wₚ₋ * ψₚ .+ net.λ₋) .* ψₜₚ

