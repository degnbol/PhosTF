#!/usr/bin/env julia
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
isdefined(Main, :FluxUtils) || include("../utilities/FluxUtils.jl")

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

struct Mdl
    nₜ::Int # number of TFs
    nₚ::Int # number of KPs
    nₒ::Int # number of non-regulators
    Wₜ # raw trainable weights for TF->V regularization. Shape nᵥ×nₜ+nₚ due to efficiency of multiplication in model. Real weights are found after applying Mₜ.
    Wₚ # raw trainable weights for KP->TFKP regularization. Shape nₜ+nₚ×nₜ+nₚ since it has to be square for the inversion. Real weights are found after applying Mₚ.
    Mₜ::Matrix{Bool} # mask to restrict which elements of Wₜ that are nonzero 
    Mₚ::Matrix{Bool} # mask to restrict which elements of Wₚ that are nonzero
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end
struct MdlS
    nₜ::Int # number of TFs
    nₚ::Int # number of KPs
    nₒ::Int # number of non-regulators
    Wₜ # raw trainable weights for TF->V regularization. Real weights are found after applying Mₜ and Sₜ.
    Wₚ # raw trainable weights for KP->TFKP regularization. Real weights are found after applying Mₚ and Sₚ.
    Mₜ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₜ that are nonzero 
    Mₚ::Matrix{Bool} # nᵥ×nᵥ mask to restrict which elements of Wₚ that are nonzero
    Sₜ # sign restriction matrix for Wₜ with 0 to indicating no sign restriction, -1 to enforce negative sign and +1 to enforce positive.
    Sₚ # sign restriction matrix for Wₚ with 0 to indicating no sign restriction, -1 to enforce negative sign and +1 to enforce positive.
    U # Matrix holding column vectors of non-KO indexes for each experiment. Has to have same shape as X.
    I # The identity matrix in B = Wₜ(I - Wₚ)⁻¹
end

"- K: number of experiments (k in 1:K)"
get_model(nₜ::Integer, nₚ::Integer, nₒ::Integer, K::Integer; Wₜ=nothing, Wₚ=nothing, Mₜ=nothing, Mₚ=nothing, Sₜ=nothing, Sₚ=nothing) = begin
    get_model(nₜ, nₚ, nₒ, I(nₜ+nₚ+nₒ, K); Wₜ=Wₜ, Wₚ=Wₚ, Mₜ=Mₜ, Mₚ=Mₚ, Sₜ=Sₜ, Sₚ=Sₚ)
end
"- J: Matrix holding column vectors of KO indexes for each experiment. Same shape as X."
get_model(nₜ::Integer, nₚ::Integer, nₒ::Integer, J::AbstractMatrix; Wₜ=nothing, Wₚ=nothing, Mₜ=nothing, Mₚ=nothing, Sₜ=nothing, Sₚ=nothing) = begin
    Wₜ !== nothing || (Wₜ = random_weight(nₜ+nₚ+nₒ, nₜ))
    Wₚ !== nothing || (Wₚ = random_weight(nₜ+nₚ, nₚ))
    Mₜ !== nothing || (Mₜ = _Mₜ(nₜ+nₚ+nₒ, nₜ))
    Mₚ !== nothing || (Mₚ = _Mₚ(nₜ, nₚ))
    Wₜ = [Wₜ zeros(nₜ+nₚ+nₒ, nₚ)]
    Wₚ = [zeros(nₜ+nₚ, nₜ) Wₚ]
    Mₜ = [Mₜ zeros(Bool, nₜ+nₚ+nₒ, nₚ)]
    Mₚ = [zeros(Bool, nₜ+nₚ, nₜ) Mₚ]
    _get_model(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ, Mₚ, J, Sₜ, Sₚ)
end
_get_model(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix, Sₜ::AbstractMatrix, Sₚ::AbstractMatrix) = MdlS(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ, Mₚ, Sₜ, Sₚ, J2U(J), eye(Mₚ))
_get_model(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix, Sₜ::Nothing, Sₚ::Nothing) = _get_model(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ, Mₚ, J)
_get_model(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ::Matrix{Bool}, Mₚ::Matrix{Bool}, J::AbstractMatrix) = Mdl(nₜ, nₚ, nₒ, Wₜ, Wₚ, Mₜ, Mₚ, J2U(J), eye(Mₚ))

"Make the struct models trainable either using both Wₜ and Wₚ or only the latter."
function make_trainable(train_Wₜ::Bool=true)
    @functor Mdl
    @functor MdlS
    
    if train_Wₜ
        @eval Flux.trainable(m::Mdl) = (m.Wₜ, m.Wₚ)
        @eval Flux.trainable(m::MdlS) = (m.Wₜ, m.Wₚ)
    else
        @eval Flux.trainable(m::Mdl) = (m.Wₚ,)
        @eval Flux.trainable(m::MdlS) = (m.Wₚ,)
    end
end


"Restrict sign of entries in W with a matrix S where 0 means no restriction, -1 means entry has to have negative sign and +1 for positive."
apply_S(W, S) = W .* (S .== 0) .+ abs.(W) .* S

"Get weights corrected for masks, etc. Note that Wₜ has shape nᵥ×nₜ+nₚ so dimensions match for internal model multiplication."
_Wₜ(m::Mdl) = m.Wₜ .* m.Mₜ
_Wₚ(m::Mdl) = m.Wₚ .* m.Mₚ
_Wₜ(m::MdlS) = apply_S(m.Wₜ, m.Sₜ) .* m.Mₜ
_Wₚ(m::MdlS) = apply_S(m.Wₚ, m.Sₚ) .* m.Mₚ
"Get the weights in like above except the shape of Wₜ is nᵥ×nₜ"
WₜWₚ(m) = _Wₜ(m)[:, 1:m.nₜ], _Wₚ(m)[:, m.nₜ+1:m.nₜ+m.nₚ]

""" "predict" logFC X given itself using the model. x̂ = B x + e, where B = Wₜ(I - Wₚ)⁻¹"""
(m::Mdl)(x) = _Wₜ(m) * ((I - _Wₚ(m)) \ x)
(m::MdlS)(x) = _Wₜ(m) * ((I - _Wₚ(m)) \ x)

_B(m) = _Wₜ(m) * inv(I - _Wₚ(m))

# B* is B calculated from Wₜ and Wₚ with absolute elements.
_Bstar(m) = abs.(_Wₜ(m)) * inv(I - abs.(_Wₚ(m)))

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
"- ks: If we are using batches, then indicate which batches are used"
SSE(m, X::Matrix, ks) = sum(abs2, (m.U[:, ks] .* (m(X) .- X)))


"Loss function to train parameters in W to result in a B that is as similar to a solution to B from LLC method (Eberhardt)."
SSE_B(m, B_LLC::Matrix) = sum(abs2, (_B(m) .- B_LLC))

SSE_T(m, X::Matrix) = begin
	n, K = size(X)
	sum(abs2, ([X2T(X) zeros(n, n-K)] - _T(m)))
end

end;
