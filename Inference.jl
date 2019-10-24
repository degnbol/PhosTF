#!/usr/bin/env julia
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end
if !isdefined(Main, :FluxUtils) include("utilities/FluxUtils.jl") end

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Statistics: mean
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
using Printf
using ..ArrayUtils: eye, shuffle_columns
include("Model.jl"); using .Model
import ..FluxUtils

function random_W(n::Integer, m::Integer)
	out = FluxUtils.random_weight(n::Int, m::Int)
	out[diagind(out)] .= 0
	out
end

"W is the param weight matrix, W′ is the masked version where untrainable entries are set to zero."
L1(X, W, W′, constants, λ::Real) = Model.sse(constants, W′, X) + λ*Model.l1(W)
L1(X, W, W′, Ex, Ey, constants, λ::Real) = Model.sse(constants, W′, X, Ex, Ey) + λ*Model.l1(W)
L1_alt(X, W, W′, constants, λ::Real) = Model.sse_alt(constants, W′, X) + λ*Model.l1(W)
L1_alt(X, W, W′, Ex, Ey, constants, λ::Real) = Model.sse_alt(constants, W′, X, Ex, Ey) + λ*Model.l1(W)

loss = [
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> L1(X, W, W′, constants, λ::Real),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + .1Model.l1(W)),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W)),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> Model.sse(cs, W′, X) + λ*(Model.l1(Model._T(Model._B(cs, W))) + .1Model.l1(W)),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> Model.sse(cs, W′, X) + λ*(Model.l1(Model._T(abs.(Model._B(cs, W)))) + .1Model.l1(W)),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> Model.sse(cs, W′, X) + λ*(Model.l1(Model._T(Model._B(cs, W))) + Model.l1(W)),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> Model.sse(cs, W′, X) + λ*(Model.l1(Model._T(abs.(Model._B(cs, W)))) + Model.l1(W)),
]

loss_alt = [
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> L1_alt(X, W, W′, constants, λ::Real),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse_alt(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W)),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> Model.sse_alt(constants, W′, X) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W) + Model.l1(W*Iₜ))
]

loss_EE = [
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> L1(X, W, W′, Ex, Ey, constants, λ::Real),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X, Ex, Ey) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W) + Model.l1(Ex) + Model.l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> Model.sse(constants, W′, X, Ex, Ey) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W) + Model.l1(W*Iₜ) + Model.l1(Ex) + Model.l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> L1_alt(X, W, W′, Ex, Ey, constants, λ::Real),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> Model.sse_alt(constants, W′, X, Ex, Ey) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W) + Model.l1(Ex) + Model.l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> Model.sse_alt(constants, W′, X, Ex, Ey) + λ*(Model.l_cas(W′, Iₜ, Iₚ) + Model.l1(W) + Model.l1(W*Iₜ) + Model.l1(Ex) + Model.l1(Ex))
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
		l = L(X)/length(W)
		mse = mean(abs.(Model.sse(constants, Model.apply_priors(W, M, S), X)))
		@printf("%.5f\t%.5f\n", l, mse)
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	collect(Model.apply_priors(W, nothing, S))
end

function infer_alt(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing, conf=0)
	n, K = size(X)
	constants = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	Ex = param(FluxUtils.random_weight(n, K))
	Ey = param(FluxUtils.random_weight(n, K))
	if M == nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	
	L(X) = loss_EE[conf](X, W, Model.apply_priors(W, M, S), Ex, Ey, constants, Model.Iₜ(n, nₜ, nₚ), Model.Iₚ(n, nₜ, nₚ), λ)
	
	function cb()
		l = L(X)/length(W)
		ex, ey = mean(abs.(Ex)), mean(abs.(Ey))
		mse = mean(abs.(Model.sse(constants, Model.apply_priors(W, M, S), X, Ex, Ey)))
		@printf("%.5f\t%.5f\t%.5f\t%.5f\n", l, mse, ex, ey)
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
