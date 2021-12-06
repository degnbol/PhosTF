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
Get masks M, S for trainable weights and/or restriction to sign of weights.
- mat: Matrix{String} or BitMatrix, or other Matrix with 0(/unrecognized)=no edge, 1==possible edge, "+"==positive edge, "-"==negative edge.
return: mask M for trainable weights, mask S for sign restriction
"""
function mat2MS(mat::AbstractMatrix)
	positives = mat .== "+"
	negatives = mat .== "-"
	S = positives - negatives
    M = (mat .== 1) .| (S .!= 0)
    if !any(M .== 0) M = nothing end
    if !any(S .!= 0) S = nothing end
    M, S
end
mat2MS(::Nothing) = nothing, nothing


function infer(X, nₜ::Integer, nₚ::Integer, out_WT::String="WT_infer.mat", out_WP::String="WP_infer.mat"; J=nothing, epochs::Integer=5000, opt="ADAMW", 
        lr::Float64=0.001, decay::Real=0, WT=nothing, WP=nothing, WT_mask=nothing, WP_mask=nothing,
        lambda_Bstar::Real=.1, lambda_absW::Real=0., reg_WT::Bool=true, train_WT::Bool=true)
    # read matrices that were given
	X = loaddlm_(X, Float64)
    J = loaddlm_(J, Float64)
    initial_Wₜ = loaddlm_(WT)
    initial_Wₚ = loaddlm_(WP)
    Mₜ, Sₜ = mat2MS(loaddlm_(WT_mask))
    Mₚ, Sₚ = mat2MS(loaddlm_(WP_mask))
    
    J === nothing || @assert(size(J) == size(X))
    
	nᵥ, K = size(X)
	nₒ = nᵥ - (nₜ + nₚ)

    mdl = Model.Mdl(nₜ, nₚ, nₒ, J === nothing ? K : J; Wₜ=initial_Wₜ, Wₚ=initial_Wₚ, Mₜ=Mₜ, Mₚ=Mₚ, Sₜ=Sₜ, Sₚ=Sₚ)
    Model.make_trainable(train_WT)
	
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

