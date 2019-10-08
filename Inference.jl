#!/usr/bin/env julia
"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
include("Model.jl")
include("utilities/ArrayUtils.jl"); using .ArrayUtils: eye, shuffle_columns

loss(X) = Model.mse(constants, W, X) + norm(W)

opt = ADAMW() # or maybe NADAM
epochs = 10000
throttle = 5  # seconds between prints

nₜ, nₚ = 4, 3
K, n = nₜ+nₚ, 2(nₜ+nₚ)
constants = Model.Constants(n, nₜ, nₚ, K)
W = param(Model.random_W(n, n))

X = randn(n,K)

function callback()
	println(loss(X))
end

println(loss(X))
data = ((X,) for _ in 1:epochs)
Flux.train!(loss, [W], data, opt, cb=Flux.throttle(callback, throttle))
println(loss(X))


end;
