#!/usr/bin/env julia
include("src/utilities/ReadWrite.jl")
include("src/utilities/CLI.jl")
include("src/utilities/General.jl")
include("src/Inference.jl")
if !isdefined(Main, :ArrayUtils) include("src/utilities/ArrayUtils.jl") end


using Fire
using LinearAlgebra
using ..ReadWrite, ..ArrayUtils, ..General
using ..Model
using ..Inference, ..CLI


"""
Get priors from files with the indicators 0=no edge, 1=possible edge, "+"=positive edge, "-"=negative edge.
Can be fed nothing values, and produces nothing values when a matrix would otherwise provide no additional information.
return: priors, priors_sign
"""
function _priors(WT_prior::Union{Matrix,Nothing}, WP_prior::Union{Matrix,Nothing}, n::Integer, nₜ::Integer, nₚ::Integer)
	if WT_prior === nothing && WP_prior === nothing return nothing, nothing end
	M, S = Model.priors(WT_prior === nothing ? n : WT_prior, WP_prior === nothing ? nₚ : WP_prior)
	if all(Model._Wₜ(M,nₜ,nₚ) .== 1) && all(Model._Wₚ(M,nₜ,nₚ) .== 1) M = nothing end
	if all(S == 0) S = nothing end
	M, S
end

"""
Infer a weight matrix from logFC data.
- WT_prior/WP_prior: optionally limit Wₜ/Wₚ if they are partially known.
- lambda: regularization constant for B*
- lambdaW: regularization constant for abs(W)
- lambdaWT: whether the WT edges are regularized. Should only be used if highly trusted WT_priors are provided.
- PKPP: fname. vector, each element is -1 or 1 indicating PP or PK. 0s ignored.
- WT/WP: previous run to continue.
- WT_reg: NOT square, weights for each element in WT for regularization cost. NaN means the weight should not be allowed, so will function as masking as well.
"""
@main function infer(X, nₜ::Integer, nₚ::Integer, ot="WT_infer.mat", op="WP_infer.mat"; epochs::Integer=5000, 
	lambda::Real=.1, lambdaW::Real=0., lambdaWT::Bool=true, WT_prior=nothing, WP_prior=nothing, PKPP=nothing, WT=nothing, WP=nothing, J=nothing,
	linex::Bool=false, trainWT::Bool=true, WT_reg=nothing)
	# empty strings is the same as providing nothing.
	WT == "" && (WT = nothing)
	WP == "" && (WP = nothing)
	WT_prior == "" && (WT_prior = nothing)
	WP_prior == "" && (WP_prior = nothing)
	WT_reg == "" && (WT_reg = nothing)
	ot, op = abspath_(ot), abspath_(op)  # weird PWD issues require abs path

	# load files
	X = loaddlm(abspath_(X), Float64)
	n = size(X,1)
	if J !== nothing
		J = loaddlm(abspath_(J), Float64)
		@assert size(J) == size(X)
	end
	WT_reg === nothing || (WT_reg = loaddlm(abspath_(WT_reg)))
	WT_prior === nothing || (WT_prior = loaddlm(abspath_(WT_prior)))
	WP_prior === nothing || (WP_prior = loaddlm(abspath_(WP_prior)))
	W = (WT === nothing || WP === nothing) ? nothing : Model._W(loaddlm(abspath_(WT)), loaddlm(abspath_(WP)))

	# use NaNs from WT_reg for masking
	if WT_reg !== nothing
		nans = isnan.(WT_reg)
		if any(nans)
			if (WT_prior === nothing) WT_prior = .!nans
			# if we have a mask given explicitly then all values that were NaN should be removed by the mask
			else @assert all(WT_prior[nans] .== 0) end
			# set them to 1 (default weight) to avoid any NaN related problems.
			WT_reg[nans] = 1
		end
	end
	M, S = _priors(WT_prior, WP_prior, n, nₜ, nₚ)
	
	if PKPP !== nothing
		PKPP = vec(loaddlm(abspath_(PKPP)))
		padding = zeros(n-length(PKPP))
		Iₚₖ = diagm([PKPP .== +1; padding])
		Iₚₚ = diagm([PKPP .== -1; padding])
	else
		Iₚₖ, Iₚₚ = nothing, nothing
	end

	W_reg = WT_reg === nothing ? nothing : [ones(n,nₚ) WT_reg]

	W = Inference.infer(X, nₜ, nₚ; epochs=epochs, λ=lambda, λW=lambdaW, λWT=lambdaWT, M=M, S=S, Iₚₖ=Iₚₖ, Iₚₚ=Iₚₚ, W=W, J=J, linex=linex, trainWT=trainWT, W_reg=W_reg)
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

