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

L1(X, W, constants, λ::Real) = Model.sse(constants, W, X) + λ*Model.l1(W)
loss = [
(X, W, constants, M, Iₜ, Iₚ, λ::Real) ->
L1(X, W, constants, λ),
(X, W, constants, M, Iₜ, Iₚ, λ::Real) ->
Model.sse(constants, W, X) + λ*(.1Model.l_cas(W.*M, Iₜ, Iₚ) + Model.l1(W, Iₜ)),
(X, W, constants, M, Iₜ, Iₚ, λ::Real) ->
Model.sse(constants, W, X) + λ*(.1sum(Model.node_cas(W.*M, Iₜ)) + Model.l1(W, Iₜ))
]

"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), conf=0)
	n, K = size(X)
	constants = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	
	L(X) = loss[conf](X, W, constants, offdiag(W), Model.Iₜ(n, nₜ, nₚ), Model.Iₚ(n, nₜ, nₚ), λ)
	
	function cb()
		println(L(X)/length(W))
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	collect(W)
end
function train!(W, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
end

end;
