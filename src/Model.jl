#!/usr/bin/env julia
isdefined(Main, :ArrayUtils) || include("utilities/ArrayUtils.jl")
isdefined(Main, :FluxUtils) || include("utilities/FluxUtils.jl")

"Core model describing the equations of the model etc."
module Model
using LinearAlgebra
using Statistics: mean
using Flux
using ..ArrayUtils: eye
import ..FluxUtils

export offdiag, random_W
export sse, sse_B, sse_T
export l1
export _B, B_star, _T

"A mask to remove diagonal of a matrix."
function offdiag(matrix)
	out = ones(size(matrix))
	out[diagind(out)] .= 0
	out
end

function random_W(n::Integer)
	out = FluxUtils.random_weight(n,n)
	out[diagind(out)] .= 0
	out
end


function default_Mₜ(nᵥ::Integer, nₜ::Integer, nₚ::Integer)::Matrix{Bool}
	mat = zeros(Bool, nᵥ, nᵥ)
	mat[:, nₚ+1:nₚ+nₜ] .= 1
	mat[diagind(mat)] .= 0
	mat
end
function default_Mₚ(nᵥ::Integer, nₜ::Integer, nₚ::Integer)::Matrix{Bool}
	mat = zeros(Bool, nᵥ, nᵥ)
	mat[1:nₜ+nₚ, 1:nₚ] .= 1
	mat[diagind(mat)] .= 0
	mat
end
_U(J)::Matrix = 1 .- J

"""
Constant terms in the equations of the inference model.
- J: Matrix holding column vectors of KO indexes for each experiment. No need to be square but has to have same shape as X.
return:
- Mₜ: Masking matrix for TF. Square.
- Mₚ: Masking matrix for KP. Square.
- U: Matrix holding column vectors of non-KO indexes for each experiment. No need to be square but has to have same shape as X.
"""
constants(nᵥ::Integer, nₜ::Integer, nₚ::Integer, J::AbstractMatrix) = (Mₜ=default_Mₜ(nᵥ,nₜ,nₚ), Mₚ=default_Mₚ(nᵥ,nₜ,nₚ), U=_U(J))
"- K: number of experiments (k in 1:K)"
constants(nᵥ::Integer, nₜ::Integer, nₚ::Integer, K::Integer=nₜ+nₚ) = constants(nᵥ,nₜ,nₚ, eye(nᵥ,K))
constants(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix) = constants(nᵥnₜnₚ(Wₜ,Wₚ)...)


nᵥnₜnₚ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix) = size(Wₜ,1), size(Wₜ,2), size(Wₚ,2)
function nₒnₜnₚ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	(nᵥ,nₜ), nₚ = size(Wₜ), size(Wₚ,2)
	nᵥ-(nₜ+nₚ),nₜ,nₚ
end
nₒnₜnₚ(WₜWₚ::Array) = nₒnₜnₚ(WₜWₚ...)
nₒnₜnₚ(WₜWₚ::Tuple) = nₒnₜnₚ(WₜWₚ...)

function _W(Wₜ, Wₚ)
	nₒ, nₜ, nₚ = nₒnₜnₚ(Wₜ, Wₚ)
	[[Wₚ; zeros(nₒ,nₚ)] Wₜ zeros(nₒ+nₜ+nₚ,nₒ)]
end

_Wₜ(W::AbstractMatrix, nₜ::Integer, nₚ::Integer) = W[:,nₚ+1:nₚ+nₜ]
_Wₚ(W::AbstractMatrix, nₜ::Integer, nₚ::Integer) = W[1:nₜ+nₚ,1:nₚ]
_Wₜ(W::AbstractMatrix, Iₜ) = W*Iₜ
_Wₚ(W::AbstractMatrix, Iₚ) = W*Iₚ
_Wₜ(W::Tuple, Iₜ) = W[1]*Iₜ
_Wₚ(W::Tuple, Iₚ) = W[2]*Iₚ
WₜWₚ(W::AbstractMatrix, nₜ::Integer, nₚ::Integer) = _Wₜ(W,nₜ,nₚ), _Wₚ(W,nₜ,nₚ)
WₜWₚ(W::Tuple, nₜ::Integer, nₚ::Integer) = _Wₜ(W[1],nₜ,nₚ), _Wₚ(W[2],nₜ,nₚ)
"Assert that untrainable areas are in fact zero."
function isW(W::AbstractMatrix, nₜ::Integer, nₚ::Integer)
	all(W[nₜ+nₚ+1:end,1:nₚ] .== 0) && 
	all(W[:,nₜ+nₚ+1:end] .== 0) && 
	all(diag(W) .== 0)
end
function isW(W::Tuple, nₜ::Integer, nₚ::Integer)
	all(W[1][:,1:nₚ] .== 0) && 
	all(W[1][:,nₜ+nₚ+1:end] .== 0) && 
	all(W[2][nₜ+nₚ+1:end,1:nₚ] .== 0) && 
	all(W[2][:,nₚ+1:end] .== 0) && 
	all(diag(W[1]) .== 0) &&
	all(diag(W[2]) .== 0)
end


Iₚ(nᵥ::Integer, nₜ::Integer, nₚ::Integer) = diagm([[1 for _ in 1:nₚ]; [0 for _ in nₚ+1:nᵥ]])
Iₜ(nᵥ::Integer, nₜ::Integer, nₚ::Integer) = diagm([[0 for _ in 1:nₚ]; [1 for _ in 1:nₜ]; [0 for _ in nₚ+nₜ+1:nᵥ]])
Iₒ(nᵥ::Integer, nₜ::Integer, nₚ::Integer) = diagm([[0 for _ in 1:nₚ+nₜ]; [1 for _ in nₚ+nₜ+1:nᵥ]])
"All except nodes ∈ KP"
I₋ₖₚ(Iₚ₊::AbstractMatrix, Iₚ₋::AbstractMatrix) = I(size(Iₚ₊,1)) - (Iₚ₊+Iₚ₋)

_B(W::AbstractMatrix, cs::NamedTuple) = (W.*cs.Mₜ) * inv(I(size(W,1)) - W.*cs.Mₚ)
"""
- W[1]: Wₜ, square
- W[2]: Wₚ, square
I has to have a known size to not produce an error that might be fixed in later release.
"""
_B(W::Tuple, cs::NamedTuple) = (W[1].*cs.Mₜ) * inv(I(size(W[2],1)) - W[2].*cs.Mₚ)
"To avoid finding inverse matrix, we can instead solve if given the x in B^-1 * x"
_B(W::AbstractMatrix, cs::NamedTuple, x) = (W.*cs.Mₜ) * ((I(size(W,1)) - W.*cs.Mₚ) \ x)
_B(W::Tuple, cs::NamedTuple, x) = (W[1].*cs.Mₜ) * ((I(size(W[2],1)) - W[2].*cs.Mₚ) \ x)
B_star(W::AbstractMatrix, cs::NamedTuple) = _B(abs.(W), cs)
B_star(W::Tuple, cs::NamedTuple) = _B([abs.(W[1]), abs.(W[2])], cs)


_T(B) = (I(size(B,1)) - (B.*offdiag(B))) \ B
_T(W::AbstractMatrix, cs::NamedTuple) = _T(_B(W, cs))
"""
Get the total effects from each node to each node as a simple linear regression coefficient.
Section 6.1 of Eberhardt report "Learning Linear Cyclic Causal Models with Latent Variables".
Note that self-loops are removed.
"""
X2T(X) = X.*offdiag(X) ./ repeat(diag(X)', size(X,1), 1)

"""
Error. Difference between prediction and truth.
- W: either Matrix with Wt and Wp, or vector with Wt and Wp matrices, or tracked versions.
"""
E(W, cs::NamedTuple, X::Matrix) = cs.U .* (_B(W,cs,X) .- X)

"""
- cs: struct containing the constants Mₜ, Mₚ, and U
- W: either Matrix with Wt and Wp, or vector with Wt and Wp matrices, or tracked versions.
- X: Matrix holding column vectors of measured (simulated) logFC values. No need to be square but has to have the same shape as J.
"""
sse(W, cs::NamedTuple, X::Matrix) = sum(E(W,cs,X) .^ 2)
"- ks: If we are using batches, then indicate which batches are used"
sse(W, cs::NamedTuple, X::Matrix, ks) = sum((cs.U[:,ks] .* (_B(W,cs,X) .- X)) .^ 2)
"Alternative SSE where both TF and KP edges onto a KO are removed instead of only TF. Reduces edges among KP."
function sse_alt(W::AbstractMatrix, cs::NamedTuple, X::Matrix)
	i = I(size(W,1))
	E = cs.U .* X .- hcat([(W.*cs.Mₜ .* uₖ) * ((i - W.*cs.Mₚ .* uₖ) \ x) for (uₖ, x) in zip(eachcol(cs.U), eachcol(X))]...)
	sum(E.^2)
end


"""
Loss function to train parameters in W to result in a B that is as similar to a solution to B from LLC method (Eberhardt).
"""
sse_B(W::AbstractMatrix, cs::NamedTuple, B_LLC::Matrix) = sum((_B(W,cs) .- B_LLC).^2)

function sse_T(W::AbstractMatrix, cs::NamedTuple, X::Matrix)
	n,K = size(X)
	sum(([X2T(X) zeros(n,n-K)] - _T(W,cs)).^2)
end

l1(W::AbstractMatrix) = norm(W,1)
l1(W::Tuple) = norm(W[2],1)


"""
Get a mask for trainable weights and/or restriction to sign of weights.
- W: matrix with 0(/unrecognized)=no edge, 1==possible edge, "+"==positive edge, "-"==negative edge.
Also works if W is a BitMatrix with true (1) and false (0)
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
priors(nᵥ::Integer, Wₚ_prior::AbstractMatrix) = priors(ones(nᵥ, size(Wₚ_prior,1)-size(Wₚ_prior,2)), Wₚ_prior)
"""
Get priors from files with the indicators 0=no edge, 1=possible edge, "+"=positive edge, "-"=negative edge.
Can be fed nothing values, and produces nothing values when a matrix would otherwise provide no additional information.
- WT_prior/WP_prior: should be either matrix with 0,1,+,- or bitmatrix.
return: priors, priors_sign
"""
function priors(WT_prior::Union{AbstractMatrix,Nothing}, WP_prior::Union{AbstractMatrix,Nothing}, nᵥ::Integer, nₜ::Integer, nₚ::Integer)
	if WT_prior === nothing && WP_prior === nothing return nothing, nothing end
	M, S = Model.priors(WT_prior === nothing ? nᵥ : WT_prior, WP_prior === nothing ? nₚ : WP_prior)
	if all(Model._Wₜ(M,nₜ,nₚ) .== 1) && all(Model._Wₚ(M,nₜ,nₚ) .== 1) M = nothing end
    any(S .!= 0) || (S = nothing)
	M, S
end



apply_priors(W, M, S::AbstractMatrix) = apply_priors(apply_priors(W, M), nothing, S)
apply_priors(W::AbstractMatrix, ::Nothing, S::AbstractMatrix) = W .* (S.==0) .+ abs.(W) .* S
apply_priors(W::Tuple, ::Nothing, S::AbstractMatrix) = [apply_priors(W[1], nothing, S), apply_priors(W[2], nothing, S)]
apply_priors(W::AbstractMatrix, M::AbstractMatrix) = W.*M
apply_priors(W::Tuple, M::AbstractMatrix) = [W[1].*M,    W[2].*M]
apply_priors(W::Tuple, M::Tuple) = [W[1].*M[1], W[2].*M[2]]
apply_priors(W, M, ::Nothing) = apply_priors(W, M)
apply_priors(W, ::Nothing, ::Nothing) = W

end;
