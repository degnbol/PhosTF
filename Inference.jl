#!/usr/bin/env julia
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end

"Flux machine learning for gradient descent of errors defined by model loss functions."
module Inference
using LinearAlgebra
using Random
using Statistics: mean
using Flux, Flux.Tracker
using Flux.Tracker: grad, update!
using ..ArrayUtils: eye, shuffle_columns
include("Model.jl"); using .Model
using Formatting

"""
Get random indexes taken from ∈ [1,K] in portions.
return: a generator where each iteration will provide a list of size "batch_size" with indexes ∈ [1,K].
An index will not appear twice. If "K" is not divisible by "batch_size" there will be randomly left out indexes that will not be returned.
"""
function batches(K, batch_size)
	ks = shuffle(1:K)
	n_batches = trunc(Int, K/batch_size)
	(ks[1+batch_size*(i-1):batch_size*i] for i ∈ 1:n_batches)
end


function train!(W, V, X::AbstractMatrix, L, cb, epochs, opt)
	data = ((X,) for epoch ∈ 1:epochs)
	params = (V === nothing) ? [W] : [W, V]
	Flux.train!(L, params, data, opt, cb=cb)
end
"Using batches. Batch size should be K when training finishes. Changing batch size not implemented."
function train!(W, V, X::AbstractMatrix, L, cb, epochs, opt, batch_size)
	data = ((X[:,ks], ks) for epoch ∈ 1:epochs for ks ∈ batches(size(X,2), batch_size))
	params = (V === nothing) ? [W] : [W, V]
	Flux.train!(L, params, data, opt, cb=cb)
end


"W is the param weight matrix, W′ is the masked version where untrainable entries are set to zero."
L1(X, W′, cs, λ::Real) = sse(cs, W′, X) + λ*l1(W)
loss(X, W, W′, cs, λ::Real) = sse(cs, W′, X) + λ*l1(_B(cs, abs.(W)))


get_V(::Nothing, ::Nothing, ::Any) = nothing
get_V(Iₚₖ::Matrix, Iₚₚ::Matrix, ::Nothing) = FluxUtils.random_weight(size(Iₚₖ,1),1) .|> abs |> param
get_V(Iₚₖ::Matrix, Iₚₚ::Matrix, W::Matrix) = sign.(sum(W*(Iₚₖ-Iₚₚ); dims=2)) |> param


"""
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
- Iₚₖ, Iₚₚ: kinase and phosphatase indicator diagonal matrices
- W: from previous training.
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ::Real=.1, throttle=5, opt=ADAMW(), 
	M=nothing, S=nothing, Iₚₖ=nothing, Iₚₚ=nothing, W=nothing)
	n, K = size(X)
	if M === nothing M = ones(n, n) end # no prior knowledge
	M[diagind(M)] .= 0  # enforce no self loops
	cs = Model.Constants(n, nₜ, nₚ, K)
	V = get_V(Iₚₖ, Iₚₚ, W)
	W !== nothing || (W = FluxUtils.random_weight(n,n))
	# mask W before tracking it, so some entries are untrainable
	W = param(W .* (cs.Mₜ .+ cs.Mₚ) .* M)
	Iₚ = V === nothing ? Model.Iₚ(n, nₜ, nₚ) : Iₚₖ + Iₚₚ
	Iₜ = Model.Iₜ(n, nₜ, nₚ)
	Iₓ = I(n) - (Iₜ+Iₚ)

	L(X) = loss(X, W, Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ), cs, λ)

	function cb()
		l = L(X)
		w′ = Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ)
		e = Model.sse(cs, w′, X)
		lt = l1(w′*Iₜ)
		lp = l1(w′*Iₚ)
		printfmtln(5, l, e, lt, lp)
	end
	
	println("loss\tsse\tLt\tLp")
	train!(W, V, X, L, Flux.throttle(cb, throttle), epochs, opt)
	Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ)
end



end;
