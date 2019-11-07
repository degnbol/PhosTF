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
L1(X, W, W′, cs, λ::Real) = sse(cs, W′, X) + λ*l1(W)
Lsim(X, W′, cs, λ::Real) = sse(cs, W′, X) + l1(l_sim(W′, λ))
Lcas(X, W′, cs, λ::Real, Iₜ, Iₚ) = sse(cs, W′, X) + l_cas(W′, Iₜ, Iₚ)
loss(X, W′, cs, λ, Iₓ) = Lsim(X, W′, cs, λ) + l1(Iₓ*W)


"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ::Real=10, throttle=5, opt=ADAMW(), M=nothing, S=nothing)
	n, K = size(X)
	if M === nothing M = ones(n, n) end # no prior knowledge
	M[diagind(M)] .= 0  # enforce no self loops
	cs = Model.Constants(n, nₜ, nₚ, K)
	W = param(random_W(n, n))
	Iₜ = Model.Iₜ(n, nₜ, nₚ)
	Iₚ = Model.Iₚ(n, nₜ, nₚ)
	Iₓ = I(size(W,1)) - (Iₜ+Iₚ)

	L(X) = loss(X, Model.apply_priors(W, M, S), cs, λ, Iₓ)
	
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

function infer_B(B_LLC::AbstractMatrix, nₚ::Integer; epochs::Integer=10000, λ::Real=.1, throttle=5, opt=ADAMW(), M=nothing, S=nothing)
	nₚₜ, K = size(B_LLC); nₜ = nₚₜ-nₚ
	cs = Model.Constants(nₚₜ, nₜ, nₚ, K)
	W = param(random_W(nₚₜ, nₚₜ))
	if M === nothing M = offdiag(W)
	else M[diagind(M)] .= 0 end
	Iₜ = Model.Iₜ(nₚₜ, nₜ, nₚ)
	Iₚ = Model.Iₚ(nₚₜ, nₜ, nₚ)

	L(B_LLC) = loss_B(B_LLC, W, Model.apply_priors(W, M, S), cs, Iₜ, Iₚ, λ)
	
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

function train!(W, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for _ in 1:epochs)
	Flux.train!(L, [W], data, opt, cb=cb)
end

end;
