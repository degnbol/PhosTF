#!/usr/bin/env julia
include("../utilities/ReadWrite.jl")
include("../utilities/CLI.jl")
include("../utilities/ArgParseUtils.jl")
include("GradientDescent.jl")
isdefined(Main, :Model) || include("Model.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
using ArgParse
using LinearAlgebra
using Flux

# autofix_names=true changes names with dashes to underscores.
argument_parser = ArgParseSettings(description="Infer a weight matrix from logFC data.", autofix_names=true)
@add_arg_table! argument_parser begin
    "X"
        help = "Filename for LogFC values in a matrix. No column or row names. Space delimiters are recommended."
        required = true
    "nₜ"
        arg_type = Int
        range_tester = x -> x > 0
        help = "Number of TFs."
        required = true
    "nₚ"
        arg_type = Int
        range_tester = x -> x > 0
        help = "Number of KPs."
        required = true
    "out_WT"
        default = "WT_infer.mat"
        help = "Outfile for inferred Wₜ adjacency matrix."
    "out_WP"
        default = "WP_infer.mat"
        help = "Outfile for inferred Wₚ adjacency matrix."
    "--J", "-J"
        help = "Filename for J, which has 1s indicating mutated genes in each experiment and zeros otherwise. Default is using the identity matrix."
    "--epochs", "-e"
        arg_type = Int
        default = 500
        range_tester = x -> x >= 0
        help = "Number of times to train using all data in X."
    "--opt", "-o"
        default = "ADAMW"
        help = "Optimizer algorithm for gradient descent."
    "--lr", "-η"
        arg_type = Float64
        default = 0.001
        range_tester = x -> x > 0
        help = "Learning rate for gradient descent."
    "--decay", "-d"
        arg_type = Float64
        default = 0.0
        range_tester = x -> x >= 0
        help = "Weight decay for ADAMW. https://fluxml.ai/Flux.jl/stable/training/optimisers/#Flux.Optimise.ADAMW"
    "--WT", "-T"
        help = "Wₜ from a previous run to continue or to use as the true Wₜ. Default is starting from random noise."
    "--WP", "-P"
        help = "Wₚ from a previous run to continue or to use as the true Wₚ. Default is starting from random noise."
    "--WT-mask", "-t"
        help = "Masking trainable weights in Wₜ, and masking sign. 0, 1, +, and - are untrainable, trainable, positive and negative edge, respectively."
    "--WP-mask", "-p"
        help = "Masking trainable weights in Wₚ, and masking sign. 0, 1, +, and - are untrainable, trainable, positive and negative edge, respectively."
    "--lambda-Bstar", "-λ"
        arg_type = Float64
        default = 0.1
        range_tester = x -> x >= 0
        help = "Regularization factor for B*."
    "--lambda-absW", "-w"
        arg_type = Float64
        default = 0.0
        range_tester = x -> x >= 0
        help = "Regularization factor for abs(W)."
    "--reg-WT"
        arg_type = Bool
        default = true
        help = "Whether the WT edges are regularized. Should only be used if a highly trusted WT-mask is provided."
    "--train-WT"
        arg_type = Bool
        default = true
        help = "Whether the WT edges are trained (predicted) at all. If set to false, --WT/-T has to be provided with fully trusted edges."
end


function parse_optimizer(name::String, η::Float64, decay)
    if name == "ADAMW"
        # just default β
        ADAMW(η, (0.9, 0.999), decay)
    else
        # assuming learning rate is first argument, which it is for most available optimizers.
        getfield(Flux, Symbol(name))(η)
    end
end

# abspath is used throughout due to a weird PWD issue on the server.
abspath_(path::AbstractString) = abspath(expanduser(path))
loaddlm_(path::Nothing) = nothing
loaddlm_(path::AbstractString) = begin
	# empty strings is the same as providing nothing.
    path != "" || return nothing
    ReadWrite.loaddlm(abspath_(path))
end
loaddlm_(path::AbstractString, T::Type) = begin
	# empty strings is the same as providing nothing.
    path != "" || return nothing
    ReadWrite.loaddlm(abspath_(path), T)
end
savedlm_(path::AbstractString, x::AbstractArray) = begin
    ReadWrite.savedlm(abspath_(path), x)
end


"""
Get a mask for trainable weights and/or restriction to sign of weights.
- W: matrix with 0(/unrecognized)=no edge, 1==possible edge, "+"==positive edge, "-"==negative edge.
Also works if W is a BitMatrix with true (1) and false (0)
return: mask for W, mask for sign
"""
function masks(W_mask::AbstractMatrix)
	positives, negatives = W_mask .== "+", W_mask .== "-"
	possible = (W_mask .== 1) .| positives .| negatives
	possible, positives - negatives
end
function masks(Wₜ_mask::AbstractMatrix, Wₚ_mask::AbstractMatrix)
	Wₜ_mask, Wₜ_mask_sign = masks(Wₜ_mask)
	Wₚ_mask, Wₚ_mask_sign = masks(Wₚ_mask)
	_W(Wₜ_mask, Wₚ_mask), _W(Wₜ_mask_sign, Wₚ_mask_sign)
end
masks(Wₜ_mask::AbstractMatrix, nₚ::Integer) = masks(Wₜ_mask, ones(size(Wₜ_mask, 2) + nₚ, nₚ))
masks(nᵥ::Integer, Wₚ_mask::AbstractMatrix) = masks(ones(nᵥ, size(Wₚ_mask, 1) - size(Wₚ_mask, 2)), Wₚ_mask)
"""
Get masks from files with the indicators 0=no edge, 1=possible edge, "+"=positive edge, "-"=negative edge.
Can be fed nothing values, and produces nothing values when a matrix would otherwise provide no additional information.
- WT_mask/WP_mask: should be either matrix with 0,1,+,- or bitmatrix.
return: masks, masks_sign
"""
function masks(WT_mask::Union{AbstractMatrix,Nothing}, WP_mask::Union{AbstractMatrix,Nothing}, nᵥ::Integer, nₜ::Integer, nₚ::Integer)
	if WT_mask === nothing && WP_mask === nothing return nothing, nothing end
	M, S = Model.masks(WT_mask === nothing ? nᵥ : WT_mask, WP_mask === nothing ? nₚ : WP_mask)
    if all(M[:, 1:nₜ] .== 1) && all(M[1:nₜ+nₚ, nₜ+1:nₜ+nₚ]) .== 1) M = nothing end
    any(S .!= 0) || (S = nothing)
	M, S
end


function infer(X, nₜ::Integer, nₚ::Integer, out_WT::String="WT_infer.mat", out_WP::String="WP_infer.mat"; J=nothing, epochs::Integer=5000, opt="ADAMW", 
        lr::Float64=0.001, decay::Real=0, WT=nothing, WP=nothing, WT_mask=nothing, WP_mask=nothing,
        lambda_Bstar::Real=.1, lambda_absW::Real=0., reg_WT::Bool=true, train_WT::Bool=true)
    # read matrices that were given
	X = loaddlm_(X, Float64)
    J = loaddlm_(J, Float64)
    initial_Wₜ = loaddlm_(WT)
    initial_Wₚ = loaddlm_(WP)
	Wₜmask = loaddlm_(WT_mask)
	Wₚmask = loaddlm_(WP_mask)
    
    J === nothing || @assert(size(J) == size(X))
    
	nᵥ, K = size(X)
	nₒ = nᵥ - (nₜ + nₚ)
    
	M, S = masks(Wₜmask, Wₚmask, nᵥ, nₜ, nₚ)
	# TODO: apply the masks
    
	mdl = Model.Mdl(nᵥ, nₜ, nₚ, J === nothing ? K : J)
	
    # TODO: train_WT should be used to set WT weights untrainable


    opt = parse_optimizer(opt, lr, decay)
    Wₜ, Wₚ = GradientDescent.train(mdl, X; epochs=epochs, opt=opt, λBstar=lambda_Bstar, λabsW=lambda_absW, reg_Wₜ=reg_WT)
	
	# assert that code is working. If we don't intend to train WT then assert that no changes has occurred.
	if !train_WT && any(initial_Wₜ .!= Wₜ)
		n_changes = sum(initial_Wₜ .!= Wₜ)
		diff = sum(abs.(initial_Wₜ - Wₜ))
		@error("There has been $n_changes changes made to Wₜ even though it was not intented to be trained on (difference=$diff).")
		savedlm_(out_WT, Wₜ)
	end
    
	savedlm_(out_WP, Wₚ)
	train_WT && savedlm_(out_WT, Wₜ)
end

# parse args if run on command line as opposed to being imported
if abspath(PROGRAM_FILE) == @__FILE__
    ArgParseUtils.main(argument_parser, infer)
end

