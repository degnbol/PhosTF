#!/usr/bin/env julia

# totally bad quickfix
try import Flux catch; end

include("../utilities/ReadWrite.jl")
include("../utilities/CLI.jl")
include("../utilities/ArgParseUtils.jl")
include("GradientDescent.jl")
isdefined(Main, :Model) || include("Model.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")

# ArgParse used instead of Fire since Fire has a weird bug where it says there's too many arguments.
using ArgParse
using LinearAlgebra
using Flux

function argument_parser()
	s = ArgParseSettings(description="Infer a weight matrix from logFC data.", autofix_names=true)
	@add_arg_table! s begin
		"X"
			help = "Filname for LogFC values in a matrix. No column or row names. Space delimiters are recommended."
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
		"--WT-prior", "-t"
			help = "Masking trainable weights in Wₜ, and masking sign. 0, 1, +, and - are untrainable, trainable, positive and negative edge."
		"--WP-prior", "-p"
			help = "Masking trainable weights in Wₚ, and masking sign. 0, 1, +, and - are untrainable, trainable, positive and negative edge."
		"--lambda", "-λ"
			arg_type = Float64
			default = 0.1
			range_tester = x -> x >= 0
			help = "Regularization factor for B*."
		"--lambdaW", "-w"
			arg_type = Float64
			default = 0.0
			range_tester = x -> x >= 0
			help = "Regularization factor for abs(W)."
		"--lambdaWT"
			arg_type = Bool
			default = true
			help = "Whether the WT edges are regularized. Should only be used if highly trusted WT_priors are provided."
		"--trainWT"
			arg_type = Bool
			default = true
			help = "Whether the WT edges are trained at all. If set to false, --WT/-T has to be provided with fully trusted edges."
		"--WT-reg"
			help = "Filename of matrix with regularization weights for Wₜ. NOT square. 
			NaN means the weight should not be allowed, so this will function as masking as well."
	end
	s
end


function parse_optimizer(name::String, η::Float64, decay)
    if name == "ADAMW" ADAMW(η, (0.9, 0.999), decay)
    else getfield(Flux, Symbol(name))(η)  # assuming learning rate is first argument, which it is for most available optimizers. 
    end
end

abspath_(path::AbstractString) = abspath(expanduser(path))

function infer(X, nₜ::Integer, nₚ::Integer, ot::String="WT_infer.mat", op::String="WP_infer.mat"; J=nothing, epochs::Integer=5000, opt="ADAMW", lr::Float64=0.001, decay::Real=0,
	WT=nothing, WP=nothing, WT_prior=nothing, WP_prior=nothing,
	lambda::Real=.1, lambdaW::Real=0., lambdaWT::Bool=true, trainWT::Bool=true, WT_reg=nothing)
	# empty strings is the same as providing nothing.
	WT == "" && (WT = nothing)
	WP == "" && (WP = nothing)
	WT_prior == "" && (WT_prior = nothing)
	WP_prior == "" && (WP_prior = nothing)
	WT_reg == "" && (WT_reg = nothing)
	ot, op = abspath_(ot), abspath_(op)  # weird PWD issues require abs path

	# load files
	X = ReadWrite.loaddlm(abspath_(X), Float64)
	nᵥ = size(X,1)
	if J !== nothing
		J = ReadWrite.loaddlm(abspath_(J), Float64)
		@assert size(J) == size(X)
	end
	WT_reg === nothing || (WT_reg = ReadWrite.loaddlm(abspath_(WT_reg)))
	WT_prior === nothing || (WT_prior = ReadWrite.loaddlm(abspath_(WT_prior)))
	WP_prior === nothing || (WP_prior = ReadWrite.loaddlm(abspath_(WP_prior)))
	W = Model.random_W(nᵥ)
	WT = WT === nothing ? Model._Wₜ(W,nₜ,nₚ) : ReadWrite.loaddlm(abspath_(WT))
	WP = WP === nothing ? Model._Wₚ(W,nₜ,nₚ) : ReadWrite.loaddlm(abspath_(WP))
	W = Model._W(WT, WP)

	# use NaNs from WT_reg for masking
	if WT_reg !== nothing
		nans = isnan.(WT_reg)
		if any(nans)
			if (WT_prior === nothing) WT_prior = .!nans
			# if we have a mask given explicitly then all values that were NaN should be removed by the mask
			else @assert all(WT_prior[nans] .== 0) end
			# set them to 1 (default weight) to avoid any NaN related problems.
			WT_reg[nans] .= 1
		end
	end
	M, S = Model.priors(WT_prior, WP_prior, nᵥ, nₜ, nₚ)
	
	nₒ = nᵥ-(nₚ+nₜ)
	W_reg = WT_reg === nothing ? nothing : [ones(nᵥ,nₚ) WT_reg ones(nᵥ,nₒ)]

    W = Inference.infer(X, nₜ, nₚ; epochs=epochs, opt=parse_optimizer(opt, lr, decay), λ=lambda, λW=lambdaW, λWT=lambdaWT, M=M, S=S, W=W, J=J, trainWT=trainWT, W_reg=W_reg)
	Model.isW(W, nₜ, nₚ) || @error("W has nonzeros in entries that should be zero.")
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	
	if !trainWT && any(WT .!= Wₜ)
		n_changes = sum(WT .!= Wₜ)
		diff = sum(abs.(WT - Wₜ))
		@error("There has been $n_changes changes made to Wₜ even though it was not intented to be trained on (difference=$diff).")
		ReadWrite.savedlm(ot, Wₜ)
	end

	ReadWrite.savedlm(op, Wₚ)
	trainWT && ReadWrite.savedlm(ot, Wₜ)
end

# parse args if run on command line as opposed to being imported
if abspath(PROGRAM_FILE) == @__FILE__
    ArgParseUtils.main(argument_parser(), infer)
end

