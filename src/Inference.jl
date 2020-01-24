#!/usr/bin/env julia
isdefined(Main, :ArrayUtils) || include("utilities/ArrayUtils.jl")

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
using Dates

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


"W is either [Wₜ, Wₚ] or W if type is not explicit."
L1(W, ::Nothing) = 0
L1(W::AbstractMatrix, λ::Real) = λ*l1(W)
"Assuming here that if we are splitting WT and WP it means we only care to train on WP, therefore, we only reg on WP."
L1(W::AbstractVector, λ::Real) = λ*l1(W[2])
"Loss from B star"
LB(W, cs::NamedTuple, λ::Real) = λ*l1(B_star(W, cs))
LB(W, cs::NamedTuple, λ::Real, ::Nothing) = LB(W, cs, λ)
LB(W, cs::NamedTuple, λ::Real, weights) = λ*l1(B_star(W, cs).*weights)
"Only measure loss on WP not WT."
LB_WP(W, cs::NamedTuple, λ::Real) = λ*l1(B_star(W, cs).*cs.Mₚ)
LB_WP(W, cs::NamedTuple, λ::Real, ::Nothing) = LB_WP(W, cs, λ)
LB_WP(W, cs::NamedTuple, λ::Real, weights) = λ*l1(B_star(W, cs).*cs.Mₚ.*weights)

get_V(::Nothing, ::Nothing, ::Any) = nothing
get_V(Iₚₖ::Matrix, Iₚₚ::Matrix, ::Nothing) = FluxUtils.random_weight(size(Iₚₖ,1),1) .|> abs |> param
get_V(Iₚₖ::Matrix, Iₚₚ::Matrix, W::Matrix) = sign.(sum(W*(Iₚₖ-Iₚₚ); dims=2)) |> param


"""
- X: logFC values with measured nodes along axis=1, and different experiment or replicate along axis=2
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
- Iₚₖ, Iₚₚ: kinase and phosphatase indicator diagonal matrices
- W: from previous training.
- J: matrix with 1 for KO and 0 for passive observed node. Shape like X.
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ::Real=.1, λW=0., λWT::Bool=true, opt=ADAMW(), 
	M=nothing, S=nothing, Iₚₖ=nothing, Iₚₚ=nothing, W=nothing, J=nothing, linex::Bool=false, trainWT::Bool=true, W_reg=nothing)
	n, K = size(X)
	M !== nothing || (M = ones(n, n)) # no prior knowledge
	M[diagind(M)] .= 0  # enforce no self loops
	cs = Model.constants(n, nₜ, nₚ, J === nothing ? K : J)
	V = get_V(Iₚₖ, Iₚₚ, W)
	W !== nothing || (W = random_W(n))
	W = trainWT ? param(W) : [W.*cs.Mₜ, param(W.*cs.Mₚ)]
	Iₚ = V === nothing ? Model.Iₚ(n, nₜ, nₚ) : Iₚₖ + Iₚₚ
	Iₜ = Model.Iₜ(n, nₜ, nₚ)
	Iₓ = I(n) - (Iₜ+Iₚ)
	λW != 0 || (λW = nothing)
	
	error_cost = linex ? Model.linex : Model.sse
	B_cost = λWT ? LB : LB_WP
	L(X) = error_cost(Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ), cs, X) + B_cost(W, cs, λ, W_reg) + L1(W, λW)

	epoch = 0
	function cb()
		l = L(X)
		W′ = Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ)
		e = Model.sse(W′, cs, X)
		lp = l1(Model._Wₚ(W′,Iₚ))
		
		try trainWT ? printfmt(5, l,e,l1(Model._Wₜ(W′,Iₜ)),lp) : printfmt(5, l,e,lp)
		catch exception
			if isa(exception, InexactError) # if the values are crazy big we get issues with showing floats
				trainWT ? print("$(l)\t$(e)\t$(l1(Model._Wₜ(W′,Iₜ)))\t$(lp)") : print("$(l)\t$(e)\t$(lp)")
			else rethrow() end
		end
		println("\t$epoch\t$(Dates.now())")
		epoch += 1
	end
	# throttle callbacks if we are doing a small example
	n > 100 || (cb = Flux.throttle(cb, 5))
	
	println(trainWT ? "loss\tsse\tLt\tLp\tepoch\ttime" : "loss\tsse\tLp\tepoch\ttime")
	cb() # epoch 0 print before we start
	Flux.train!(L, get_params(W, V), ((X,) for _ ∈ 1:epochs), opt; cb=cb)
	Model.apply_priors(W, V, M, S, Iₚₖ, Iₚₚ)
end

get_params(W::AbstractMatrix) = [W]
get_params(W::AbstractVector) = [W[2]]  # assuming only the WP is a param
get_params(W, ::Nothing) = get_params(W)
get_params(W, V) = [get_params(W)..., V]



end;