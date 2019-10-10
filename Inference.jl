#!/usr/bin/env julia
include("utilities/ArrayUtils.jl")
include("utilities/FluxUtils.jl")

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
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

mse(X, W, constants) = Model.mse(constants, FluxUtils.zerodiag(W), X)
L1(X, W, constants, λ::Real) = mse(X, W, constants) + λ*norm(W, 1)/length(W)
function L_cas(X, W, constants, λ::Real)
	mse(X, W, constants) + 0 # TODO
end

"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs=10000, λ=.1, throttle=5, opt=ADAMW())
	n, K = size(X)
	constants = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	train!(W, X, constants, Flux.throttle(()->println(L1(X, W, constants, λ)), throttle), epochs, λ, opt)
	collect(W)
end
function train!(W, X::AbstractMatrix, constants, cb, epochs, λ, opt)
	L(X) = L1(X, W, constants, λ)
	println(L(X))
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
	println(L(X))
	nothing
end

end;
