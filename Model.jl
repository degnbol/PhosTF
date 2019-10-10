#!/usr/bin/env julia
include("utilities/ArrayUtils.jl")
include("utilities/FluxUtils.jl")

"Core model describing the equations of the model etc."
module Model
using LinearAlgebra
using Statistics: mean
using Flux
using ..ArrayUtils: eye

"Constant terms in the equations of the inference model."
struct Constants
	# Masking matrix for TF. Square.
	Mₜ::Matrix
	# Masking matrix for PK. Square.
	Mₚ::Matrix
	# Matrix holding column vectors of non-KO indexes for each experiment. No need to be square but has to have same shape as X.
	U::Matrix
	"""
	- J: Matrix holding column vectors of KO indexes for each experiment. No need to be square but has to have same shape as X.
	"""
	Constants(n::Int, nₜ::Int, nₚ::Int, J) = new(default_Mₜ(n, nₜ, nₚ), default_Mₚ(n, nₜ, nₚ), default_U(J))
	"- K: number of experiments (k in 1:K)"
	Constants(n::Int, nₜ::Int, nₚ::Int, K::Int) = Constants(n, nₜ, nₚ, eye(n,K))
	function default_Mₜ(n::Int, nₜ::Int, nₚ::Int)
		mat = zeros(Bool, n, n)
		mat[:, nₚ+1:nₚ+nₜ] .= 1
		mat[diagind(mat)] .= 0
		mat
	end
	function default_Mₚ(n::Int, nₜ::Int, nₚ::Int)
		mat = zeros(Bool, n, n)
		mat[1:nₜ+nₚ, 1:nₚ] .= 1
		mat[diagind(mat)] .= 0
		mat
	end
	default_U(J) = 1 .- J
end

nnₜnₚ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix) = size(Wₜ,1), size(Wₜ,2), size(Wₚ,2)
function nₓnₜnₚ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	(n,nₜ), nₚ = size(Wₜ), size(Wₚ,2)
	n-(nₜ+nₚ),nₜ,nₚ
end
nₓnₜnₚ(WₜWₚ::Array) = nₓnₜnₚ(WₜWₚ...)
nₓnₜnₚ(WₜWₚ::Tuple) = nₓnₜnₚ(WₜWₚ...)

function _W(Wₜ, Wₚ)
	nₓ, nₜ, nₚ = nₓnₜnₚ(Wₜ, Wₚ)
	[[Wₚ; zeros(nₓ,nₚ)] Wₜ zeros(nₓ+nₜ+nₚ,nₓ)]
end

WₜWₚ(W, nₜ, nₚ) = W[:,nₚ+1:nₚ+nₜ], W[1:nₜ+nₚ,1:nₚ]

Iₚ(n::Integer, nₜ::Integer, nₚ::Integer) = diagm([[1 for _ in 1:nₚ]; [0 for _ in nₚ+1:n]])
Iₜ(n::Integer, nₜ::Integer, nₚ::Integer) = diagm([[0 for _ in 1:nₚ]; [1 for _ in 1:nₚ]; [0 for _ in nₚ+nₜ+1:n]])

"I has to have a known size to not produce an error that might be fixed in later release."
_B(cs::Constants, W::AbstractMatrix) = W.*cs.Mₜ * (I(size(W,1)) - W.*cs.Mₚ)^-1

"""
The cummulative effects thorugh cascades.
Explicitly setting size of I so Flux can handle it.
"""
function node_cas(W::AbstractMatrix, Iₜ::AbstractMatrix)
	n = size(W,1)
	W′ = abs.(W')
	(I(n) - W′) \ (Iₜ * W′ * ones(n))
end
"Loss function for PK/PP cascades."
l_cas(W::AbstractMatrix, Iₜ::AbstractMatrix, Iₚ::AbstractMatrix) = sum((I(size(W,1)) - Iₚ*abs.(W')) \ vec(sum(Iₜ,dims=2)))
l1(W::AbstractMatrix, Iₜ::AbstractMatrix) = norm(W * Iₜ, 1)

"""
- cs: struct containing the constants Mₜ, Mₚ, and U
- W: Trainable parameters. Square matrix.
- X: Matrix holding column vectors of measured (simulated) logFC values. No need to be square but has to have the same shape as J.
"""
function sse(cs::Constants, W::AbstractMatrix, X::Matrix)
	Ue = (I - _B(cs,W)) * X .* cs.U
	sum(Ue.^2)
end

l1(W::AbstractMatrix) = norm(W,1)

"For iterating the model."
step(X::Matrix, W, cs::Constants, C::Matrix) = _B(cs,W) * X .* cs.U .+ C
step(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = (_B(cs,W) * X + E) .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix) = X .= _B(cs,W) * X .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = X .= (_B(cs,W) * X + E) .* cs.U .+ C


end;