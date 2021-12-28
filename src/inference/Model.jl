#!/usr/bin/env julia
SRC = readchomp(`git root`) * "/src/"
isdefined(Main, :ArrayUtils) || include(SRC * "utilities/ArrayUtils.jl")
isdefined(Main, :FluxUtils) || include(SRC * "utilities/FluxUtils.jl")

"Core model describing the equations of the model etc."
module Model
using LinearAlgebra
using Statistics: mean
using Flux
using Flux: @functor
using ..ArrayUtils
import ..FluxUtils: random_weight

export SSE, SSE_B, SSE_T
export _B, _Bstar, _T

"A mask to remove diagonal of a matrix."
function offdiag(matrix)
	out = ones(size(matrix))
	out[diagind(out)] .= 0
	out
end

# I was using Matrix{Bool} instead of BitMatrix since it is Int8 instead of bit storage. It might be better for the multiplications that these mats are used for.
"Mask matrices for Wₜ and Wₚ to enforce no self-edges."
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

"Smaller helper to convert Wₜ to Wₜ□ (square) and similarly for M, S."
function adjₜ2□(Wₜ::AbstractMatrix)
    nᵥ, nₜ = size(Wₜ)
    Wₜ□ = zeros(eltype(Wₜ), nᵥ, nᵥ)
    Wₜ□[:, 1:nₜ] .= Wₜ
    Wₜ□
end
function adjₚ2□(Wₚ::AbstractMatrix, nᵥ::Integer)
    nₜₚ, nₚ = size(Wₚ)
    nₜ = nₜₚ - nₚ
    Wₚ□ = zeros(eltype(Wₚ), nᵥ, nᵥ)
    Wₚ□[1:nₜₚ, nₜ+1:nₜₚ] .= Wₚ
    Wₚ□
end


struct Mdl
    nₜ::Int # number of TFs
    nₚ::Int # number of KPs
    nₒ::Int # number of non-regulators
    Wₜ□ # raw trainable weights for TF->V regularization. Real weights are found after applying Mₜ. Shape nᵥ×nᵥ.
    Wₚ□ # raw trainable weights for KP->TFKP regularization. Real weights are found after applying Mₚ. Shape nᵥ×nᵥ.
    Mₜ□::Matrix{Bool} # mask to restrict which elements of Wₜ that are nonzero 
    Mₚ□::Matrix{Bool} # mask to restrict which elements of Wₚ that are nonzero
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end
# try out a more efficient approach with smaller parts of W stored.
struct MdlE
    nᵥ::Int # number of vertices
    nₜ::Int # number of TFs
    nₚ::Int # number of KPs
    nₒ::Int # number of non-regulators
    Wₜ # size nᵥ×nₜ
    Wₚ # size nₜ+nₚ×nₚ
    Mₜ::Matrix{Bool} # mask to restrict which elements of Wₜ that are nonzero 
    Mₚ::Matrix{Bool} # mask to restrict which elements of Wₚ that are nonzero
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end
struct MdlS
    nₜ::Int # number of TFs
    nₚ::Int # number of KPs
    nₒ::Int # number of non-regulators
    Wₜ□ # raw trainable weights for TF->V regularization. Real weights are found after applying Mₜ and Sₜ. Shape nᵥ×nᵥ.
    Wₚ□ # raw trainable weights for KP->TFKP regularization. Real weights are found after applying Mₚ and Sₚ. Shape nᵥ×nᵥ.
    Mₜ□::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₜ that are nonzero 
    Mₚ□::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₚ that are nonzero
    Sₜ□ # sign restriction matrix for Wₜ with 0 to indicating no sign restriction, -1 to enforce negative sign and +1 to enforce positive.
    Sₚ□ # sign restriction matrix for Wₚ with 0 to indicating no sign restriction, -1 to enforce negative sign and +1 to enforce positive.
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end

"- K: number of experiments (k in 1:K)"
get_model(nₜ::Integer, nₚ::Integer, nₒ::Integer, K::Integer; Wₜ=nothing, Wₚ=nothing, Mₜ=nothing, Mₚ=nothing, Sₜ=nothing, Sₚ=nothing) = begin
    get_model(nₜ, nₚ, nₒ, I(nₜ+nₚ+nₒ, K); Wₜ=Wₜ, Wₚ=Wₚ, Mₜ=Mₜ, Mₚ=Mₚ, Sₜ=Sₜ, Sₚ=Sₚ)
end
"- J: Matrix holding column vectors of KO indexes for each experiment. Same shape as X."
get_model(nₜ::Integer, nₚ::Integer, nₒ::Integer, J::AbstractMatrix; Wₜ=nothing, Wₚ=nothing, Mₜ=nothing, Mₚ=nothing, Sₜ=nothing, Sₚ=nothing) = begin
    nᵥ = nₜ+nₚ+nₒ
    Wₜ !== nothing || (Wₜ = random_weight(nᵥ, nₜ))
    Wₚ !== nothing || (Wₚ = random_weight(nₜ+nₚ, nₚ))
    Mₜ !== nothing || (Mₜ = _Mₜ(nᵥ, nₜ))
    Mₚ !== nothing || (Mₚ = _Mₚ(nₜ, nₚ))
    # _get_model(nₜ, nₚ, nₒ, adjₜ2□(Wₜ), adjₚ2□(Wₚ, nᵥ), adjₜ2□(Mₜ), adjₚ2□(Mₚ, nᵥ), J, Sₜ, Sₚ)
    MdlE(nᵥ, nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ, Mₚ, J2U(J), eye(Wₚ))
end
_get_model(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□::Matrix{Bool}, Mₚ□::Matrix{Bool}, J::AbstractMatrix, Sₜ::AbstractMatrix, Sₚ::AbstractMatrix) = MdlS(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□, Mₚ□, adjₜ2□(Sₜ), adjₚ2□(Sₚ, nᵥ), J2U(J), eye(Mₚ□))
_get_model(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□::Matrix{Bool}, Mₚ□::Matrix{Bool}, J::AbstractMatrix, Sₜ::Nothing, Sₚ::Nothing) = _get_model(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□, Mₚ□, J)
_get_model(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□::Matrix{Bool}, Mₚ□::Matrix{Bool}, J::AbstractMatrix) = Mdl(nₜ, nₚ, nₒ, Wₜ□, Wₚ□, Mₜ□, Mₚ□, J2U(J), eye(Mₚ□))

@functor Mdl
@functor MdlE
@functor MdlS
    
Flux.trainable(m::Mdl) = (m.Wₜ□, m.Wₚ□)
Flux.trainable(m::MdlE) = (m.Wₜ, m.Wₚ)
Flux.trainable(m::MdlS) = (m.Wₜ□, m.Wₚ□)

"Make the struct models trainable only using Wₚ."
function untrainable_Wₜ()
    @eval Flux.trainable(m::Mdl) = (m.Wₚ□,)
    @eval Flux.trainable(m::MdlE) = (m.Wₚ,)
    @eval Flux.trainable(m::MdlS) = (m.Wₚ□,)
end

"Restrict sign of entries in W with a matrix S where 0 means no restriction, -1 means entry has to have negative sign and +1 for positive."
apply_S(W, S) = W .* (S .== 0) .+ abs.(W) .* S

"Get square weights corrected for masks, etc."
Wₜ□(m::Mdl) = m.Wₜ□ .* m.Mₜ□
Wₚ□(m::Mdl) = m.Wₚ□ .* m.Mₚ□
# we cannot differentiate copyto (.= and .+=) but we can concat.
Wₜ□(m::MdlE) = [m.Wₜ .* m.Mₜ zeros(m.nᵥ, m.nₚ+m.nₒ)]
Wₚ□(m::MdlE) = [zeros(m.nᵥ, m.nₜ) [m.Wₚ .* m.Mₚ; zeros(m.nₒ, m.nₚ)] zeros(m.nᵥ, m.nₒ)]
Wₜ□(m::MdlS) = apply_S(m.Wₜ□, m.Sₜ) .* m.Mₜ□
Wₚ□(m::MdlS) = apply_S(m.Wₚ□, m.Sₚ) .* m.Mₚ□
"Get the weights in like above except the shape of Wₜ is nᵥ×nₜ"
WₜWₚ(m) = Wₜ□(m)[:, 1:m.nₜ], Wₚ□(m)[1:m.nₜ+m.nₚ, m.nₜ+1:m.nₜ+m.nₚ]

""" "predict" logFC X given itself using the model. x̂ = B x + e, where B = Wₜ(I - Wₚ)⁻¹"""
(m::Mdl)(x) = Wₜ□(m) * ((I - Wₚ□(m)) \ x)
(m::MdlE)(x) = Wₜ□(m) * ((I - Wₚ□(m)) \ x)
(m::MdlS)(x) = Wₜ□(m) * ((I - Wₚ□(m)) \ x)

_B(m) = Wₜ□(m) * inv(I - Wₚ□(m))

# B* is B calculated from Wₜ and Wₚ with absolute elements.
_Bstar(m) = abs.(Wₜ□(m)) * inv(I - abs.(Wₚ□(m)))

"""
- Wₜ: nonsquare weight matrix for TF -> V with dimensions nᵥ×nₜ
- Wₚ: nonsquare weight matrix for KP -> TFKP with dimensions nₜ+nₚ×nₚ
"""
function nₜnₚnₒ(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	nᵥ, nₜ = size(Wₜ)
	nₚ = size(Wₚ, 2)
	nₜ, nₚ, nᵥ - (nₜ + nₚ)
end

"Total effects as defined in LLC."
_T(B::AbstractMatrix) = (I(size(B, 1)) - (B .* offdiag(B))) \ B
_T(m) = _B(m) |> _T
"""
Get the total effects from each node to each node as a simple linear regression coefficient.
Section 6.1 of Eberhardt report "Learning Linear Cyclic Causal Models with Latent Variables".
Note that self-loops are removed.
"""
X2T(X) = X .* offdiag(X) ./ repeat(diag(X)', size(X, 1), 1)

"""
Error/latents. Difference between prediction and truth. xₖ = UₖBxₖ + Uₖeₖ + cₖ
- m: model struct
- X: Matrix holding column vectors of measured or simulated logFC values. No need to be square but has to have the same shape as J.
"""
E(m, X::Matrix) = m.U .* (m(X) .- X)

SSE(m, X::Matrix) = sum(abs2, E(m, X))
MSE(m, X::Matrix) = mean(abs2, E(m, X))
"- ks: If we are using batches, then indicate which batches are used"
SSE(m, X::Matrix, ks) = sum(abs2, (m.U[:, ks] .* (m(X) .- X)))


"Loss function to train parameters in W to result in a B that is as similar to a solution to B from LLC method (Eberhardt)."
SSE_B(m, B_LLC::Matrix) = sum(abs2, (_B(m) .- B_LLC))

SSE_T(m, X::Matrix) = begin
	n, K = size(X)
	sum(abs2, ([X2T(X) zeros(n, n-K)] - _T(m)))
end

L1(x::AbstractMatrix) = sum(abs, x)
L1(m::Union{Mdl,MdlE,MdlS}) = sum(L1(W) for W in Flux.trainable(m))

end;
