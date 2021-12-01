#!/usr/bin/env julia
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
isdefined(Main, :FluxUtils) || include("../utilities/FluxUtils.jl")

"Core model describing the equations of the model etc."
module Model
using LinearAlgebra
using Statistics: mean
using Flux
using ..ArrayUtils
import ..FluxUtils: random_weight

export SSE, SSE_B, SSE_T
export _B, Bstar, _T

"A mask to remove diagonal of a matrix."
function offdiag(matrix)
	out = ones(size(matrix))
	out[diagind(out)] .= 0
	out
end

# I was using Matrix{Bool} instead of BitMatrix since it is Int8 instead of bit storage. It might be better for the multiplications that these mats are used for.
"Mask matrices for Wₜ and Wₚ to enforce no self-edges. For the square (□) model they also restrict regions of W to Wₜ and Wₚ."
function _Mₜ□(nᵥ::Integer, nₜ::Integer)::Matrix{Bool}
	mat = zeros(Bool, nᵥ, nᵥ)
	mat[:, 1:nₜ] .= 1
	mat[diagind(mat)] .= 0
	mat
end
function _Mₚ□(nᵥ::Integer, nₜ::Integer, nₚ::Integer)::Matrix{Bool}
	mat = zeros(Bool, nᵥ, nᵥ)
	mat[1:nₜ+nₚ, nₜ+1:nₜ+nₚ] .= 1
	mat[diagind(mat)] .= 0
	mat
end
function _Mₜ(nᵥ::Integer, nₜ::Integer)::Matrix{Bool}
	mat = ones(Bool, nᵥ, nₜ)
	mat[diagind(mat)] .= 0
	mat
end
function _Mₚ(nₜ::Integer, nₚ::Integer)::Matrix{Bool}
	mat = ones(Bool, nₜ+nₚ, nₚ)
	mat[diagind(mat, -nₜ)] .= 0
	mat
end

J2U(J)::Matrix = 1 .- J

"Convenience function used in model constructors to make random weights from a mask of same size."
masked_W(M) = random_weight(size(M)) .* M
"Convenience function used in model constructors to make identity matrices."
LinearAlgebra.I(mat::AbstractMatrix) = I(size(mat)...)

"A square nᵥ×nᵥ weight matrix containing both Wₜ and Wₚ."
function _W(Wₜ, Wₚ)
    @assert size(Wₜ, 1) != size(Wₜ, 2)
    @assert size(Wₚ, 1) != size(Wₚ, 2)
	nₜ, nₚ, nₒ = nₜnₚnₒ(Wₜ, Wₚ)
	[[Wₚ; zeros(nₒ, nₚ)] Wₜ zeros(nₜ+nₚ+nₒ, nₒ)]
end


struct Mdl
    Wₜ # potentially trainable weights for TF->V regularization
    Wₚ # trainable weights for KP->TFKP regularization
    Mₜ::Matrix{Bool} # mask to restrict which elements of Wₜ that are trainable 
    Mₚ::Matrix{Bool} # mask to restrict which elements of Wₚ that are trainable
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end
struct Mdl□
    W # nᵥ×nᵥ square weight matrix with source nodes in order TF, KP and target nodes in order TF, KP, O.
    nₜ::Integer # number of TFs
    nₚ::Integer # number of KPs
    Mₜ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₜ that are trainable 
    Mₚ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₚ that are trainable
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end
struct Mdl□S
    W # nᵥ×nᵥ square weight matrix with source nodes in order TF, KP and target nodes in order TF, KP, O.
    nₜ::Integer # number of TFs
    nₚ::Integer # number of KPs
    Mₜ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₜ that are trainable 
    Mₚ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₚ that are trainable
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
    S # sign restriction matrix with 0 to indicating no sign restriction, -1 to enforce negative sign and +1 to enforce positive.
end

"""
- J: Matrix holding column vectors of KO indexes for each experiment. Same shape as X.
"""
Mdl(Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix) = Mdl(masked_weight(Mₜ), masked_weight(Mₚ), Mₜ, Mₚ, J2U(J), I(Mₚ))
Mdl□(nₜ::Integer, nₚ::Integer, Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix) = Mdl□(masked_weight(Mₜ), masked_weight(Mₚ), Mₜ, Mₚ, J2U(J), I(Mₚ))
Mdl□S(nₜ::Integer, nₚ::Integer, Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix, S) = Mdl□(masked_weight(Mₜ), masked_weight(Mₚ), Mₜ, Mₚ, J2U(J), I(Mₚ), S)
Mdl(nᵥ::Integer, nₜ::Integer, nₚ::Integer, J::AbstractMatrix) = Mdl(_Mₜ(nᵥ, nₜ), _Mₚ(nₜ, nₚ), J)
Mdl□(nᵥ::Integer, nₜ::Integer, nₚ::Integer, J::AbstractMatrix) = Mdl□(nₜ, nₚ, _Mₜ□(nᵥ, nₜ), _Mₚ□(nᵥ, nₜ, nₚ), J)
Mdl□S(nᵥ::Integer, nₜ::Integer, nₚ::Integer, J::AbstractMatrix, S) = Mdl□S(nₜ, nₚ, _Mₜ□(nᵥ, nₜ), _Mₚ□(nᵥ, nₜ, nₚ), J, S)
"- K: number of experiments (k in 1:K)"
Mdl(nᵥ::Integer, nₜ::Integer, nₚ::Integer, K::Integer=nₜ+nₚ) = Mdl(nᵥ, nₜ, nₚ, I(nᵥ, K))
Mdl□(nᵥ::Integer, nₜ::Integer, nₚ::Integer, K::Integer=nₜ+nₚ) = Mdl□(nᵥ, nₜ, nₚ, I(nᵥ, K))
Mdl□S(nᵥ::Integer, nₜ::Integer, nₚ::Integer, K::Integer=nₜ+nₚ) = Mdl□S(nᵥ, nₜ, nₚ, I(nᵥ, K), S)
"""
- Wₜ: TF -> V weights. Should not be nᵥ×nᵥ.
- Wₚ: KP -> TFKP weights. Should not be nᵥ×nᵥ.
"""
Mdl(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix, J) = begin
    @assert size(Wₜ, 1) != size(Wₜ, 2)
    @assert size(Wₚ, 1) != size(Wₚ, 2)
    (nᵥ, nₜ), nₚ = size(Wₜ), size(Wₚ, 2)
    Mdl(Wₜ, Wₚ, _Mₜ(nᵥ, nₜ), _Mₚ(nₜ, nₚ), J2U(J))
end
Mdl□(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix, J) = begin
    @assert size(Wₜ, 1) != size(Wₜ, 2)
    @assert size(Wₚ, 1) != size(Wₚ, 2)
    (nᵥ, nₜ), nₚ = size(Wₜ), size(Wₚ, 2)
    Mdl□(_W(Wₜ, Wₚ), _Mₜ□(nᵥ, nₜ), _Mₚ□(nₜ, nₚ), J2U(J))
end
Mdl□S(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix, J, S) = begin
    @assert size(Wₜ, 1) != size(Wₜ, 2)
    @assert size(Wₚ, 1) != size(Wₚ, 2)
    (nᵥ, nₜ), nₚ = size(Wₜ), size(Wₚ, 2)
    Mdl□(_W(Wₜ, Wₚ), _Mₜ□(nᵥ, nₜ), _Mₚ□(nₜ, nₚ), J2U(J), S)
end


"Restrict sign of entries in W with a matrix S where 0 means no restriction, -1 means entry has to have negative sign and +1 for positive."
apply_S(W, S) = W .* (S .== 0) .+ abs.(W) .* S

"NOTE: trying I instead of mdl.I first to see if it works. Otherwise use mdl.I."
# "predict" logFC X given itself using the model. x̂ = B x + e, where B = Wₜ(I - Wₚ)⁻¹
(mdl::Mdl)(x) = (mdl.Wₜ .* mdl.Mₜ) * ((I - mdl.Wₚ .* mdl.Mₚ) \ x)
(mdl::Mdl□)(x) = (mdl.W .* mdl.Mₜ) * ((I - mdl.W .* mdl.Mₚ) \ x)
(mdl::Mdl□S)(x) = begin
    W = apply_S(mdl.W, mdl.S)
    (W .* mdl.Mₜ) * ((I - W .* mdl.Mₚ) \ x)
end
_B(mdl::Mdl) = (mdl.Wₜ .* mdl.Mₜ) * inv(I - mdl.Wₚ .* mdl.Mₚ)
_B(mdl::Mdl□) = (mdl.W .* mdl.Mₜ) * inv(I - mdl.W .* mdl.Mₚ)
_B(mdl::Mdl□S) = begin
    W = apply_S(mdl.W, mdl.S)
    (W .* mdl.Mₜ) * inv(I - W .* mdl.Mₚ)
end
# B* is B calculated from Wₜ and Wₚ with absolute elements.
Bstar(mdl::Mdl) = abs.(mdl.Wₜ .* mdl.Mₜ) * inv(I - abs.(mdl.Wₚ .* mdl.Mₚ))
Bstar(mdl::Mdl□) = abs.(mdl.W .* mdl.Mₜ) * inv(I - abs.(mdl.W .* mdl.Mₚ))
Bstar(mdl::Mdl□S) = begin
    W = apply_S(mdl.W, mdl.S)
    abs.(W .* mdl.Mₜ) * inv(I - abs.(W .* mdl.Mₚ))
end

"""
Used after training to get the resulting weights.
Applying the masks shouldn't make a difference really since they were applied to the initial noise,
however those weights are still trainable so it doesn't hurt to enforce.
"""
WₜWₚ(mdl::Mdl) = mld.Wₜ .* mdl.Mₜ, mdl.Wₚ .* mdl.Mₚ
WₜWₚ(mdl::Mdl□) = (mdl.W .* mdl.Mₜ)[:, 1:nₜ], (mdl.W .* mdl.Mₚ)[1:nₜ+nₚ, nₜ+1:nₜ+nₚ]
WₜWₚ(mdl::Mdl□S) = begin
    W = apply_S(mdl.W, mdl.S)
    W[:, 1:nₜ], W[1:nₜ+nₚ, nₜ+1:nₜ+nₚ]
end

"""
- Wₜ: nonsquare weight matrix for TF -> V
- Wₚ: nonsquare weight matrix for KP -> TFKP
"""
function nₜnₚnₒ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	(nᵥ, nₜ), nₚ = size(Wₜ), size(Wₚ, 2)
	nₜ, nₚ, nᵥ - (nₜ+nₚ)
end

"Total effects as defined in LLC."
_T(B::AbstractMatrix) = (I(size(B, 1)) - (B .* offdiag(B))) \ B
_T(mdl) = _B(mdl) |> _T
"""
Get the total effects from each node to each node as a simple linear regression coefficient.
Section 6.1 of Eberhardt report "Learning Linear Cyclic Causal Models with Latent Variables".
Note that self-loops are removed.
"""
X2T(X) = X .* offdiag(X) ./ repeat(diag(X)', size(X, 1), 1)

"""
Error/latents. Difference between prediction and truth. xₖ = UₖBxₖ + Uₖeₖ + cₖ
- mdl: model struct
- X: Matrix holding column vectors of measured (simulated) logFC values. No need to be square but has to have the same shape as J.
"""
E(mdl, X::Matrix) = mdl.U .* (mdl(X) .- X)

SSE(mdl, X::Matrix) = sum(abs2, E(mdl, X))
"- ks: If we are using batches, then indicate which batches are used"
SSE(mdl, X::Matrix, ks) = sum(abs2, (mdl.U[:, ks] .* (mdl(X) .- X)))


"Loss function to train parameters in W to result in a B that is as similar to a solution to B from LLC method (Eberhardt)."
SSE_B(mdl, B_LLC::Matrix) = sum(abs2, (_B(mdl) .- B_LLC))

SSE_T(mdl, X::Matrix) = begin
	n, K = size(X)
	sum(abs2, ([X2T(X) zeros(n, n-K)] - _T(mdl)))
end

end;
