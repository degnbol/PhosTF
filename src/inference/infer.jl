#!/usr/bin/env julia
# totally bad quickfix
#= try import Flux catch; end =#

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
    "ot"
        default = "WT_infer.mat"
        help = "Outfile for inferred Wₜ adjacency matrix."
    "op"
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
    "--lambda", "-λ"
        arg_type = Float64
        default = 0.1
        range_tester = x -> x >= 0
        help = "Regularization factor for B*."
    "--lambda-W", "-w"
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



function infer(X, nₜ::Integer, nₚ::Integer, ot::String="WT_infer.mat", op::String="WP_infer.mat"; J=nothing, epochs::Integer=5000, opt="ADAMW", 
        lr::Float64=0.001, decay::Real=0, WT=nothing, WP=nothing, WT_mask=nothing, WP_mask=nothing,
        lambda::Real=.1, lambda_W::Real=0., reg_WT::Bool=true, train_WT::Bool=true)
    # read matrices that were given
	X = loaddlm_(X, Float64)
    J = loaddlm_(J, Float64)
    Wₜ = loaddlm_(WT)
    Wₚ = loaddlm_(WP)
	Wₜmask = loaddlm_(WT_mask)
	Wₚmask = loaddlm_(WP_mask)
    
    J === nothing || @assert(size(J) == size(X))
    
	nᵥ = size(X, 1)
	nₒ = nᵥ - (nₚ + nₜ)
    
    # set to random noise if not provided
    Wₜ !== nothing || (Wₜ = Model._Wₜ(Model.random_W(nᵥ), nₜ, nₚ))
    Wₚ !== nothing || (Wₚ = Model._Wₚ(Model.random_W(nᵥ), nₜ, nₚ))
	W = Model._W(Wₜ, Wₚ)
    
	M, S = Model.masks(Wₜmask, Wₚmask, nᵥ, nₜ, nₚ)
    
    opt = parse_optimizer(opt, lr, decay)
    W = GradientDescent.infer(X, nₜ, nₚ; epochs=epochs, opt=opt, λ=lambda, λW=lambda_W, λWT=lambda_WT, M=M, S=S, W=W, J=J, train_WT=train_WT)
	# assert that code is working.
	Model.isW(W, nₜ, nₚ) || @error("W has nonzeros in entries that should be zero.")
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	
	# assert that code is working. If we don't intend to train WT then assert that no changes has occurred.
	if !train_WT && any(WT .!= Wₜ)
		n_changes = sum(WT .!= Wₜ)
		diff = sum(abs.(WT - Wₜ))
		@error("There has been $n_changes changes made to Wₜ even though it was not intented to be trained on (difference=$diff).")
		savedlm_(ot, Wₜ)
	end
    
	savedlm_(op, Wₚ)
	train_WT && savedlm_(ot, Wₜ)
end

# parse args if run on command line as opposed to being imported
if abspath(PROGRAM_FILE) == @__FILE__
    ArgParseUtils.main(argument_parser, infer)
end

