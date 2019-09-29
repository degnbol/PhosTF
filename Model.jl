#!/usr/bin/env julia
module Model
using LinearAlgebra
using Statistics: mean
using Flux
include("utilities/ArrayUtils.jl"); using .ArrayUtils: eye
include("utilities/FluxUtils.jl")

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
		mat[:, 1:nₜ] .= 1
		mat[diagind(mat)] .= 0
		mat
	end
	function default_Mₚ(n::Int, nₜ::Int, nₚ::Int)
		mat = zeros(Bool, n, n)
		mat[1:nₜ+nₚ, nₜ+1:nₜ+nₚ] .= 1
		mat[diagind(mat)] .= 0
		mat
	end
	default_U(J) = 1 .- J
end

function _W(Wₜ, Wₚ)
	(n, nₜ), nₚ = size(Wₜ), size(Wₚ,2)
	nₒ = n-(nₜ+nₚ)
	[Wₜ [Wₚ; zeros(nₒ,nₚ)] zeros(n,nₒ)]
end

random_W(n::Int, m::Int) = FluxUtils.zerodiag(FluxUtils.random_weight(n::Int, m::Int))

"I has to have a known size to not produce an error that might be fixed in later release."
_B(cs::Constants, W::T where T<:AbstractMatrix) = W.*cs.Mₜ * (I(size(W,1)) - W.*cs.Mₚ)^-1

"""
- cs: struct containing the constants Mₜ, Mₚ, and U
- W: Trainable parameters. Square matrix.
- X: Matrix holding column vectors of measured (simulated) logFC values. No need to be square but has to have the same shape as J.
"""
function mse(cs::Constants, W::T where T<:AbstractMatrix, X::Matrix)
	Ue = (I - _B(cs,W)) * X .* cs.U
	mean(Ue.^2)
end

step(X::Matrix, W, cs::Constants, C::Matrix) = _B(cs,W) * X .* cs.U .+ C
step(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = (_B(cs,W) * X + E) .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix) = X .= _B(cs,W) * X .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = X .= (_B(cs,W) * X + E) .* cs.U .+ C


end;