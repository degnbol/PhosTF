#!/usr/bin/env julia
include("utilities/ArrayUtils.jl")
include("utilities/FluxUtils.jl")

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Statistics: mean
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
using ..ArrayUtils: eye, shuffle_columns
include("Model.jl"); using .Model
import ..FluxUtils

function random_W(n::Integer, m::Integer)
	out = FluxUtils.random_weight(n::Int, m::Int)
	out[diagind(out)] .= 0
	out
end

"A mask to remove diagonal of a matrix."
function offdiag(matrix)
	out = ones(size(matrix))
	out[diagind(out)] .= 0
	out
end

"W is the param weight matrix, W′ is the masked version where untrainable entries are set to zero."
L1(X, W, W′, constants, λ::Real) = Model.sse(constants, W′, X) + λ*Model.l1(W)
loss = [
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> L1(X, W, W′, constants, λ::Real),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W*Iₜ)),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W*Iₜ) + Model.l1(W)),
]


"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing, conf=0)
	n, K = size(X)
	constants = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	if M == nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	
	L(X) = loss[conf](X, W, Model.apply_priors(W, M, S), constants, Model.Iₜ(n, nₜ, nₚ), Model.Iₚ(n, nₜ, nₚ), λ)
	
	function cb()
		println(L(X)/length(W))
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	collect(Model.apply_priors(W, nothing, S))
end
function train!(W, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
end

end;
