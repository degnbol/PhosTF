#!/usr/bin/env julia
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
isdefined(Main, :ReadWrite) || include("../utilities/ReadWrite.jl")
isdefined(Main, :Model) || include("Model.jl")

"Flux machine learning for gradient descent of errors defined by model loss functions."
module GradientDescent
using LinearAlgebra
using Random
using Statistics: mean
using Formatting
using Dates
using Flux
using ..ArrayUtils: eye, shuffle_columns
using ..ReadWrite
using ..Model

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


"W is either (Wₜ, Wₚ) or W if type is not explicit."
L1(W, ::Nothing) = 0
L1(W::AbstractMatrix, λ::Real) = λ*l1(W)
"Assuming here that if we are splitting WT and WP it means we only care to train on WP, therefore, we only reg on WP."
L1(W::Tuple, λ::Real) = λ*l1(W[2])
"Loss from B star"
LB(W, cs::NamedTuple, λ::Real) = λ*l1(B_star(W, cs))
LB(W, cs::NamedTuple, λ::Real, ::Nothing) = LB(W, cs, λ)
LB(W, cs::NamedTuple, λ::Real, weights) = λ*l1(B_star(W, cs).*weights)
"Only measure loss on WP not WT."
LB_WP(W, cs::NamedTuple, λ::Real) = λ*l1(B_star(W, cs).*cs.Mₚ)
LB_WP(W, cs::NamedTuple, λ::Real, ::Nothing) = LB_WP(W, cs, λ)
LB_WP(W, cs::NamedTuple, λ::Real, weights) = λ*l1(B_star(W, cs).*cs.Mₚ.*weights)


"""
- X: logFC values with measured nodes along axis=1, and different experiment or replicate along axis=2
- throttle: seconds between prints
- opt: ADAMW or maybe NADAM
- W: from previous training.
- J: matrix with 1 for KO and 0 for passive observed node. Shape like X.
- λWT: should Wₜ be regularized on?
- trainWT: should Wₜ be trained or only Wₚ? if both are trained it is done together using W, otherwise we use W=[Wₜ,param(Wₚ)]
- save_every: e.g. 10 to save every tenth epoch. Use zero to not save intermediates. Intermediates are saved to W{T,P}.mat.tmp in PWD.
"""
function infer(X::AbstractMatrix, nₜ::Integer, nₚ::Integer; epochs::Integer=10000, λ::Real=.1, λW=0., λWT::Bool=true, opt=ADAMW(), 
	M=nothing, S=nothing, W=nothing, J=nothing, trainWT::Bool=true, W_reg=nothing, save_every::Integer=1)
	@assert W !== nothing
	nᵥ, K = size(X)
	cs = Model.constants(nᵥ, nₜ, nₚ, J === nothing ? K : J)
	M === nothing && (M = ones(nᵥ,nᵥ)) # no prior knowledge
	M = trainWT ? M .* (cs.Mₜ .+ cs.Mₚ) : (M .* cs.Mₜ, M .* cs.Mₚ) # enforce masks
	W = trainWT ? W : (W .* cs.Mₜ, W .* cs.Mₚ)
	Iₚ = Model.Iₚ(nᵥ, nₜ, nₚ)
	Iₜ = Model.Iₜ(nᵥ, nₜ, nₚ)
	λW == 0 && (λW = nothing)
	
	B_cost = λWT ? LB : LB_WP
	function L(X)
		W′ = Model.apply_priors(W, M, S)
		Model.sse(W′, cs, X) + B_cost(W′, cs, λ, W_reg) + L1(W′, λW)
	end
	
	epoch = 0
	function cb()
		l = L(X)
		W′ = Model.apply_priors(W, M, S)
		e = Model.sse(W′, cs, X)
		lp = l1(Model._Wₚ(W′,Iₚ))
		
		try trainWT ? printfmt(5, l,e,l1(Model._Wₜ(W′,Iₜ)),lp) : printfmt(5, l,e,lp)
		catch exception
			if isa(exception, InexactError) # if the values are crazy big we get issues with showing floats
				trainWT ? print("$l\t$e\t$(l1(Model._Wₜ(W′,Iₜ)))\t$lp") : print("$l\t$e\t$lp")
			else rethrow() end
		end
		println("\t$epoch\t$(Dates.now())")
		epoch += 1
		
		if save_every > 0 && epoch % save_every == 0
			W′ = Model.apply_priors(W, M, S)
			Model.isW(W′, nₜ, nₚ) || @error("W has nonzeros in entries that should be zero")
			Wₜ, Wₚ = Model.WₜWₚ(W′, nₜ, nₚ)
			trainWT && savedlm("WT.tmp.mat", Wₜ)
			savedlm("WP.tmp.mat", Wₚ)
		end
	end
	# throttle callbacks if we are doing a small example.
    # BUG: epoch is wrong when throttled, it counts number of times cb is called.
	nᵥ > 100 || (cb = Flux.throttle(cb, 5))
	
	println(trainWT ? "loss\tsse\tLt\tLp\tepoch\ttime" : "loss\tsse\tLp\tepoch\ttime")
	cb() # epoch 0 print before we start
	Flux.train!(L, get_params(W), ((X,) for _ ∈ 1:epochs), opt; cb=cb)
	Model.apply_priors(W, M, S)
end

get_params(W::AbstractMatrix) = [W]
get_params(W::Tuple) = [W[2]]  # assuming only the WP is a param



end;
