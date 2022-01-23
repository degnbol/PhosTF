#!/usr/bin/env julia
"Flux machine learning for gradient descent of errors defined by model loss functions."
module GradientDescent
using LinearAlgebra
using Random
using Statistics: mean
using Printf
using Dates
using Flux
Main.@use "utilities/ReadWrite"
Main.@use "inference/Model"
using ..Model: L1, MSE

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
        return X -> SSE(mdl, X) + λBstar * L1(_Bstar(mdl)) + λabsW * L1(mdl)
    else
        return X -> SSE(mdl, X) + λBstar * L1(_Bstar(mdl) .* mdl.Mₚ) + λabsW * L1(mdl)
    end
end
get_loss_func(mdl, λBstar::Float64, λabsW::Nothing, reg_Wₜ::Bool) = begin
	if reg_Wₜ
        return X -> SSE(mdl, X) + λBstar * L1(_Bstar(mdl))
    else
        return X -> SSE(mdl, X) + λBstar * L1(_Bstar(mdl) .* mdl.Mₚ)
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
- X: logFC values with measured nodes along axis=1 sorted in order TF, KP, O, and different experiment or replicate along axis=2
- epochs: number of epochs to train for
- λBstar: reg factor for B*
- λabsW: reg for abs(W)
- opt: ADAMW or maybe NADAM
- W: from previous training.
- J: matrix with 1 for KO and 0 for passive observed node. Shape like X.
- reg_Wₜ: should Wₜ be regularized on?
- save_times: e.g. 10 to save 10 intermediary results 10 times under gradient descent. Intermediates are saved to W{T,P}.mat.tmp in PWD.
"""
function train(mdl, X::AbstractMatrix, log::IO=stdout; epochs::Integer=10000, λBstar::Real=.1, λabsW::Real=0., opt=ADAMW(), reg_Wₜ::Bool=true, save_times::Integer=0)
    # Flux.trainable(mdl) will be set to either (Wₚ,) or (Wₜ, Wₚ) in Model, so we use it to see if we intent to train Wₜ.
    train_Wₜ = length(Flux.trainable(mdl)) == 2
    loss = get_loss_func(mdl, λBstar, λabsW, reg_Wₜ)
    save_every = save_times == 0 ? 0 : Int(epochs / save_times)

	function cb(epoch::Integer)
        line = @sprintf "%.3f\t%.3f\t%.3f" loss(X) MSE(mdl, X) L1(Model.Wₚ□(mdl))
        if train_Wₜ line *= @sprintf "\t%.3f" L1(Model.Wₜ□(mdl)) end
        line *= "\t$epoch\t$(Dates.now())"
		println(log, line)
		
		if save_times > 0 && epoch % save_every == 0
            Wₜ, Wₚ = Model.WₜWₚ(mdl)
			train_Wₜ && savedlm("WT.tmp.mat", Wₜ)
			savedlm("WP.tmp.mat", Wₚ)
		end
	end
	# throttle callbacks if we are doing a small example.
	_cb = size(X, 1) > 100 ? cb : Flux.throttle(cb, 5)
	
    println(log, train_Wₜ ? "loss\tMSE\tL1(Wp)\tL1(Wt)\tepoch\ttime" : "loss\tMSE\tL1(Wp)\tepoch\ttime")
	_cb(0) # epoch 0 print before we start
    for epoch in 1:epochs
        Flux.train!(loss, params(mdl), ((X,),), opt)
        _cb(epoch)
    end
    cb(epochs) # last epoch print after finish
    Model.WₜWₚ(mdl)
end
# if a log file name is provided
function train(mdl, X::AbstractMatrix, log::AbstractString; kwargs...)
    open(log, "w") do fh
        return train(mdl, X, fh; kwargs...)
    end
end

end;
