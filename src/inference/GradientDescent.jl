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
function batches(K::Integer, batch_size::Integer)
	ks = shuffle(1:K)
	n_batches = trunc(Int, K/batch_size)
	(ks[1+batch_size*(i-1):batch_size*i] for i ∈ 1:n_batches)
end


# If any of the λs are zero then avoid adding that part to the loss, e.g. using ::Nothing
get_loss_func(mdl, λBstar::Float64, λabsW::Float64, reg_Wₜ::Bool) = begin
    if λBstar == 0 return get_loss_func(mdl, nothing, λabsW, reg_Wₜ) end
    if λabsW == 0 return get_loss_func(mdl, λabsW, nothing, reg_Wₜ) end
	if reg_Wₜ
        return X -> SSE(mdl, X) + λBstar * L1(Bstar(mdl)) + λabsW * L1(mdl)
    else
        return X -> SSE(mdl, X) + λBstar * L1(Bstar(mdl) .* mdl.Mₚ) + λabsW * L1(mdl)
    end
end
get_loss_func(mdl, λBstar::Float64, λabsW::Nothing, reg_Wₜ::Bool) = begin
	if reg_Wₜ
        return X -> SSE(mdl, X) + λBstar * L1(Bstar(mdl))
    else
        return X -> SSE(mdl, X) + λBstar * L1(Bstar(mdl) .* mdl.Mₚ)
    end
end
# TODO: should reg_Wₜ have effect on the λabsW term? Also in the Float64, Float64 func
get_loss_func(mdl, λBstar::Nothing, λabsW::Float64, reg_Wₜ::Bool) = begin
    if reg_Wₜ
        return X -> SSE(mdl, X) + λabsW * L1(mdl)
    else
        return X -> SSE(mdl, X) + λabsW * L1(mdl)
    end
end


"""
- mdl: Model
- X: logFC values with measured nodes along axis=1, and different experiment or replicate along axis=2
- epochs: number of epochs to train for
- λBstar: reg factor for B*
- λabsW: reg for abs(W)
- opt: ADAMW or maybe NADAM
- W: from previous training.
- J: matrix with 1 for KO and 0 for passive observed node. Shape like X.
- reg_Wₜ: should Wₜ be regularized on?
- train_Wₜ: should Wₜ be trained or only Wₚ? if both are trained it is done together using W, otherwise we use W=[Wₜ,param(Wₚ)]
- save_every: e.g. 10 to save every tenth epoch. Use zero to not save intermediates. Intermediates are saved to W{T,P}.mat.tmp in PWD.
"""
function train(mdl, X::AbstractMatrix; epochs::Integer=10000, λBstar::Real=.1, λabsW::Real=0., opt=ADAMW(), reg_Wₜ::Bool=true, save_every::Integer=1)
    loss = get_loss_func(mdl, λBstar, λabsW, reg_Wₜ)

	epoch = 0
	function cb()
		l = loss(X)
		e = SSE(mdl, X)
		lp = L1(Model._Wₚ(mdl))
		
		try train_Wₜ ? printfmt(5, l, e, L1(Model._Wₜ(mdl)), lp) : printfmt(5, l, e, lp)
		catch exception
			if isa(exception, InexactError) # if the values are crazy big we get issues with showing floats
				train_Wₜ ? print("$l\t$e\t$(L1(Model._Wₜ(W)))\t$lp") : print("$l\t$e\t$lp")
			else rethrow() end
		end
		println("\t$epoch\t$(Dates.now())")
		epoch += 1
		
		if save_every > 0 && epoch % save_every == 0
			Wₜ, Wₚ = Model.WₜWₚ(mdl)
			train_WT && savedlm("WT.tmp.mat", Wₜ)
			savedlm("WP.tmp.mat", Wₚ)
		end
	end
	# throttle callbacks if we are doing a small example.
    # BUG: epoch is wrong when throttled, it counts number of times cb is called.
	nᵥ = size(X, 1)
	nᵥ > 100 || (cb = Flux.throttle(cb, 5))
	
	println(train_WT ? "loss\tSSE\tLt\tLp\tepoch\ttime" : "loss\tSSE\tLp\tepoch\ttime")
	cb() # epoch 0 print before we start
	Flux.train!(loss, params(mdl), ((X,) for _ ∈ 1:epochs), opt; cb=cb)
    WₜWₚ(mdl)
end

end;
