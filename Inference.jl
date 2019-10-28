#!/usr/bin/env julia
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end
if !isdefined(Main, :FluxUtils) include("utilities/FluxUtils.jl") end

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Statistics: mean
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
using ..ArrayUtils: eye, shuffle_columns
include("Model.jl"); using .Model
using Formatting
import ..FluxUtils

function random_W(n::Integer, m::Integer)
	out = FluxUtils.random_weight(n::Int, m::Int)
	out[diagind(out)] .= 0
	out
end

"W is the param weight matrix, W′ is the masked version where untrainable entries are set to zero."
L1(X, W, W′, cs, λ::Real) = sse(cs, W′, X) + λ*l1(W)
L1(X, W, W′, Ex, Ey, cs, λ::Real) = sse(cs, W′, X, Ex, Ey) + λ*l1(W)
L1_alt(X, W, W′, cs, λ::Real) = sse_alt(cs, W′, X) + λ*l1(W)
L1_alt(X, W, W′, Ex, Ey, cs, λ::Real) = sse_alt(cs, W′, X, Ex, Ey) + λ*l1(W)

loss = [
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse(cs, W′, X) + .1l1(W*Iₜ) + .01l1(W*Iₚ),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse(cs, W′, X) + .1l_cas(W′, Iₜ, Iₚ) + .1l1(W*Iₜ) + .01l1(W*Iₚ),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse(cs, W′, X) + .1l1(_T(_B(cs, W′))) + .1l1(W*Iₜ) + .01l1(W*Iₚ),
(X, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse(cs, W′, X) + .1(_B(cs, W′) .|> abs |> _T |> l1) + .1l1(W*Iₜ) + .01l1(W*Iₚ),
]
loss_B = [
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B(cs, W′, B_LLC) + .2l1(W*Iₜ) + .01l1(W*Iₚ),
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B(cs, W′, B_LLC) + .1l_cas(W′, Iₜ, Iₚ) + .2l1(W*Iₜ) + .01l1(W*Iₚ),
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B(cs, W′, B_LLC) + .1(_B(cs, W′) .|> abs |> _T |> l1) + .2l1(W*Iₜ) + .01l1(W*Iₚ),
]
loss_B_alt = [
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B_alt(cs, W′, B_LLC) + .1l1(W),
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B_alt(cs, W′, B_LLC) + l_cas(W′, Iₜ, Iₚ) + .1l1(W),
(B_LLC, W, W′, cs, Iₜ, Iₚ, λ::Real) -> sse_B_alt(cs, W′, B_LLC) + l1(_T(_B_alt(cs, W′))) + .1l1(W),
]
loss_alt = [
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> L1_alt(X, W, W′, constants, λ::Real),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> sse_alt(constants, W′, X) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W)),
(X, W, W′, constants, Iₜ, Iₚ, λ::Real) -> sse_alt(constants, W′, X) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W) + l1(W*Iₜ))
]
loss_EE = [
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> L1(X, W, W′, Ex, Ey, constants, λ::Real),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> sse(constants, W′, X, Ex, Ey) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W) + l1(Ex) + l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> sse(constants, W′, X, Ex, Ey) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W) + l1(W*Iₜ) + l1(Ex) + l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> L1_alt(X, W, W′, Ex, Ey, constants, λ::Real),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> sse_alt(constants, W′, X, Ex, Ey) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W) + l1(Ex) + l1(Ex)),
(X, W, W′, Ex, Ey, constants, Iₜ, Iₚ, λ::Real) -> sse_alt(constants, W′, X, Ex, Ey) + λ*(l_cas(W′, Iₜ, Iₚ) + l1(W) + l1(W*Iₜ) + l1(Ex) + l1(Ex))
]


"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing, conf=0)
	n, K = size(X)
	cs = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	if M === nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	Iₜ = Model.Iₜ(n, nₜ, nₚ)
	Iₚ = Model.Iₚ(n, nₜ, nₚ)

	L(X) = loss[conf](X, W, Model.apply_priors(W, M, S), cs, Iₜ, Iₚ, λ)
	
	function cb()
		l = L(X)
		w′ = Model.apply_priors(W, M, S)
		e = Model.sse(cs, w′, X)
		lt = l1(W*Iₜ)
		lp = l1(W*Iₚ)
		cas = l_cas(w′, Iₜ, Iₚ)
		T = l1(_T(_B(cs, w′)))
		T_abs = l1(_T(abs.(_B(cs, w′))))
		printfmtln(5, l, e, lt, lp, cas, T, T_abs)
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	
	println("loss\tsse\tLt\tLp\tcas\tT\tT_abs")
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	Model.apply_priors(W, nothing, S)
end

function infer_B(B_LLC::AbstractMatrix, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing, conf=1)
	nₚₜ, K = size(B_LLC); nₜ = nₚₜ-nₚ
	cs = Model.Constants(nₚₜ, nₜ, nₚ, K)
	W = param(random_W(nₚₜ, nₚₜ))
	if M === nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	Iₜ = Model.Iₜ(nₚₜ, nₜ, nₚ)
	Iₚ = Model.Iₚ(nₚₜ, nₜ, nₚ)

	L(B_LLC) = loss_B[conf](B_LLC, W, Model.apply_priors(W, M, S), cs, Iₜ, Iₚ, λ)
	
	function cb()
		l = L(B_LLC)
		w′ = Model.apply_priors(W, M, S)
		sse = sse_B(cs, w′, B_LLC)
		lt = l1(W*Iₜ)
		lp = l1(W*Iₚ)
		cas = l_cas(w′, Iₜ, Iₚ)
		T = l1(_T(_B(cs, w′)))
		T_abs = l1(_T(abs.(_B(cs, w′))))
		printfmtln(5, l, sse, lt, lp, cas, T, T_abs)
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end

	println("loss\tsse\tLt\tLp\tcas\tT\tT_abs")
	train!(W, B_LLC, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	Model.apply_priors(W, nothing, S)
end

function infer_EE(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing, conf=0)
	n, K = size(X)
	constants = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	Ex = param(FluxUtils.random_weight(n, K))
	Ey = param(FluxUtils.random_weight(n, K))
	if M === nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	
	L(X) = loss_EE[conf](X, W, Model.apply_priors(W, M, S), Ex, Ey, constants, Model.Iₜ(n, nₜ, nₚ), Model.Iₚ(n, nₜ, nₚ), λ)
	
	function cb()
		l = L(X)/length(W)
		mse = mean(abs.(sse(constants, Model.apply_priors(W, M, S), X, Ex, Ey)))
		ex, ey = mean(abs.(Ex)), mean(abs.(Ey))
		printfmtln(5, l, mse, ex, ey)
		d = diag(W); if any(abs.(d) .> .001) Flux.error("Nonzero diagonal") end
	end
	
	train!(W, X, L, Flux.throttle(cb, throttle), epochs, opt)
	# not using M just to make it clear if values supposed to be zero for some reason are not (sanity check)
	Model.apply_priors(W, nothing, S)
end


function train!(W, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
end

end;
