#!/usr/bin/env julia
include("utilities/ArrayUtils.jl")
include("utilities/FluxUtils.jl")

"Core model describing the equations of the model etc."
module Model
using LinearAlgebra
using Statistics: mean
using Flux
using ..ArrayUtils: eye
import ..FluxUtils

export offdiag, random_W
export sse, sse2, sse_B, sse_T
export l1
export _B, _T

"A mask to remove diagonal of a matrix."
function offdiag(matrix)
	out = ones(size(matrix))
	out[diagind(out)] .= 0
	out
end


"Constant terms in the equations of the inference model."
struct Constants
	# Masking matrix for TF. Square.
	Mₜ::Matrix{Bool}
	# Masking matrix for PK. Square.
	Mₚ::Matrix{Bool}
	# Matrix holding column vectors of non-KO indexes for each experiment. No need to be square but has to have same shape as X.
	U::Matrix
	"""
	- J: Matrix holding column vectors of KO indexes for each experiment. No need to be square but has to have same shape as X.
	"""
	Constants(n::Int, nₜ::Int, nₚ::Int, J) = new(default_Mₜ(n, nₜ, nₚ), default_Mₚ(n, nₜ, nₚ), _U(J))
	"- K: number of experiments (k in 1:K)"
	Constants(n::Int, nₜ::Int, nₚ::Int, K::Int=nₜ+nₚ) = Constants(n, nₜ, nₚ, eye(n,K))
end

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
_U(J) = 1 .- J

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

_Wₜ(W, nₜ::Integer, nₚ::Integer) = W[:,nₚ+1:nₚ+nₜ]
_Wₚ(W, nₜ::Integer, nₚ::Integer) = W[1:nₜ+nₚ,1:nₚ]
"""
wₚᵢⱼ := vᵢ⋅(|wₚₖᵢⱼ| - |wₚₚᵢⱼ|)
As vector notation:
Wₚ := V (Wₚₖ' - Wₚₚ')
We use rowwise multiplication with .* here instead of diagonal matrix multiplication so make sure V is a 2D column vector.
"""
_Wₚ(W, V, Iₚₖ::Matrix, Iₚₚ::Matrix) = V .* (abs.(W)*(Iₚₖ-Iₚₚ))
WₜWₚ(W, nₜ::Integer, nₚ::Integer) = _Wₜ(W,nₜ,nₚ), _Wₚ(W,nₜ,nₚ)


Iₚ(n::Integer, nₜ::Integer, nₚ::Integer) = diagm([[1 for _ in 1:nₚ]; [0 for _ in nₚ+1:n]])
Iₜ(n::Integer, nₜ::Integer, nₚ::Integer) = diagm([[0 for _ in 1:nₚ]; [1 for _ in 1:nₜ]; [0 for _ in nₚ+nₜ+1:n]])

"I has to have a known size to not produce an error that might be fixed in later release."
_B(cs::Constants, W::AbstractMatrix) = (W.*cs.Mₜ) * inv(I(size(W,1)) - W.*cs.Mₚ)
"To avoid finding inverse matrix, we can instead solve if given the x in B^-1 * x"
_B(cs::Constants, W::AbstractMatrix, x) = (W.*cs.Mₜ) * ((I(size(W,1)) - W.*cs.Mₚ) \ x)

_T(B) = (I(size(B,1)) - (B.*offdiag(B))) \ B
_T(cs::Constants, W::AbstractMatrix) = _T(_B(cs, W))
"""
Get the total effects from each node to each node as a simple linear regression coefficient.
Section 6.1 of Eberhardt report "Learning Linear Cyclic Causal Models with Latent Variables".
Note that self-loops are removed.
"""
X2T(X) = X.*offdiag(X) ./ repeat(diag(X)', size(X,1), 1)


"""
- cs: struct containing the constants Mₜ, Mₚ, and U
- W: Trainable parameters. Square matrix.
- X: Matrix holding column vectors of measured (simulated) logFC values. No need to be square but has to have the same shape as J.
"""
function sse(cs::Constants, W::AbstractMatrix, X::Matrix)
	E = X .- _B(cs,W,X)
	sum((cs.U .* E) .^ 2)
end
"- ks: If we are using batches, then indicate which batches are used"
function sse(cs::Constants, W::AbstractMatrix, X::Matrix, ks)
	E = X .- _B(cs,W,X)
	sum((cs.U[:,ks] .* E) .^ 2)
end
"Alternative SSE where both TF and PK/PP edges onto a KO are removed instead of only TF. Reduces edges among PK/PP."
function sse_alt(cs::Constants, W::AbstractMatrix, X::Matrix)
	i = I(size(W,1))
	E = cs.U .* X .- hcat([(W.*cs.Mₜ .* uₖ) * ((i - W.*cs.Mₚ .* uₖ) \ x) for (uₖ, x) in zip(eachcol(cs.U), eachcol(X))]...)
	sum(E.^2)
end

"""
Loss function to train parameters in W to result in a B that is as similar to a solution to B from LLC method (Eberhardt).
"""
sse_B(cs::Constants, W::AbstractMatrix, B_LLC::Matrix) = sum((_B(cs,W) .- B_LLC).^2)

function sse_T(cs::Constants, W::AbstractMatrix, X::Matrix)
	n,K = size(X)
	sum(([X2T(X) zeros(n,n-K)] - _T(cs, W)).^2)
end

l1(W::AbstractMatrix) = norm(W,1)


"""
Get a mask for trainable weights and/or restriction to sign of weights.
- W: matrix with 0(/unrecognized)=no edge, 1=possible edge, "+"=positive edge, "-"=negative edge.
return: mask for W, mask for sign
"""
function priors(W_prior::AbstractMatrix)
	positives, negatives = W_prior .== "+", W_prior .== "-"
	possible = (W_prior .== 1) .| positives .| negatives
	possible, positives - negatives
end
function priors(Wₜ_prior::AbstractMatrix, Wₚ_prior::AbstractMatrix)
	Wₜ_prior, Wₜ_prior_sign = priors(Wₜ_prior)
	Wₚ_prior, Wₚ_prior_sign = priors(Wₚ_prior)
	_W(Wₜ_prior, Wₚ_prior), _W(Wₜ_prior_sign, Wₚ_prior_sign)
end
priors(Wₜ_prior::AbstractMatrix, nₚ::Integer) = priors(Wₜ_prior, ones(size(Wₜ_prior,2)+nₚ, nₚ))
priors(n::Integer, Wₚ_prior::AbstractMatrix) = priors(ones(n, size(Wₚ_prior,1)-size(Wₚ_prior,2)), Wₚ_prior)

apply_priors(W, M, S) = apply_priors(W .* M, nothing, S)
apply_priors(W, ::Nothing, S) = W .* (S .== 0) .+ abs.(W) .* S
apply_priors(W, M, ::Nothing) = W .* M
apply_priors(W, ::Nothing, ::Nothing) = W
"If we know which nodes are ∈ PK and ∈ PP, then use that information."
apply_priors(W, V::AbstractArray, M, S, Iₚₖ::Matrix, Iₚₚ::Matrix) = apply_priors(W*(I(size(W,1))-(Iₚₖ+Iₚₚ)) + _Wₚ(W,V,Iₚₖ,Iₚₚ), M, S)
apply_priors(W, ::Nothing, M, S, ::Any, ::Any) = apply_priors(W, M, S)

end;