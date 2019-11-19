#!/usr/bin/env julia
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Statistics: mean
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
using ..ArrayUtils: eye, shuffle_columns
include("Model.jl"); using .Model
using Formatting

"W is the param weight matrix, W′ is the masked version where untrainable entries are set to zero."
L1(X, W′, cs, λ::Real) = sse(cs, W′, X) + λ*l1(W′)
loss(X, W′, cs, λ_W::Real, λ_B::Real) = sse(cs, W′, X) + λ*l1(W′) + λ*l1(_B(cs, abs.(W′)))

"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ_W::Real=.1, λ_B::Real=.1, throttle=5, opt=ADAMW(), 
	M=nothing, S=nothing, Iₚₖ=nothing, Iₚₚ=nothing)
	n, K = size(X)
	if M === nothing M = ones(n, n) end # no prior knowledge
	M[diagind(M)] .= 0  # enforce no self loops
	cs = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	Iₚ = (Iₚₖ === nothing || Iₚₚ === nothing) ? Model.Iₚ(n, nₜ, nₚ) : Iₚₖ + Iₚₚ
	Iₜ = Model.Iₜ(n, nₜ, nₚ)
	Iₓ = I(n) - (Iₜ+Iₚ)

	L(X) = loss(X, Model.apply_priors(W, M, S, Iₚₖ, Iₚₚ), cs, λ_W, λ_B)

	function cb()
		l = L(X)
		w′ = Model.apply_priors(W, M, S)
		e = Model.sse(cs, w′, X)
		lt = l1(W*Iₜ)
		lp = l1(W*Iₚ)
		printfmtln(5, l, e, lt, lp)
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	
	println("loss\tsse\tLt\tLp")
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	Model.apply_priors(W, nothing, S, Iₚₖ, Iₚₚ)
end

function train!(W, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
end

end;
