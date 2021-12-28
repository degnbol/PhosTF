#!/usr/bin/env julia
SRC = readchomp(`git root`) * "/src/"
include(SRC * "utilities/ReadWrite.jl")
include(SRC * "utilities/CLI.jl")
include(SRC * "utilities/ArgParseUtils.jl")
include(SRC * "inference/GradientDescent.jl")
isdefined(Main, :Model) || include(SRC * "inference/Model.jl")
isdefined(Main, :ArrayUtils) || include(SRC * "utilities/ArrayUtils.jl")
isdefined(Main, :DataFrameUtils) || include(SRC * "utilities/DataFrameUtils.jl")
using ArgParse
using LinearAlgebra
using Flux
using DataFrames


# autofix_names=true changes names with dashes to underscores.
argument_parser = ArgParseSettings(description="Infer a weight matrix from logFC data.", autofix_names=true)
@add_arg_table! argument_parser begin
    "logFC"
        help = "Filename for LogFC values in a matrix. Column names are mutated gene(s) in an experiment, first column has row names which are measured genes."
        required = true
    "TF"
        help = "Filename with a gene name per line that are (potential) TF. Alternatively supply a regex pattern to match against mutated genes. If nothing is provided, then all mutated genes are considered potential."
    "KP"
        help = "Filename with a gene name per line that are (potential) KP. Alternatively supply a regex pattern to match against mutated genes. If nothing is provided, then all mutated genes are considered potential."
    "out_WT"
        default = "WT_infer.mat"
        help = "Outfile for inferred Wₜ adjacency matrix."
    "out_WP"
        default = "WP_infer.mat"
        help = "Outfile for inferred Wₚ adjacency matrix."
    "--log"
        default = stdout
        help = "File to write log to. Default is stdout."
    "--mut-sep", "-s"
        arg_type = Char
        help = "Character that separate gene names given in a column name when multiple genes are mutated in an experiment. By default all column names are assumed to be a single gene."
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
readlines_(path::AbstractString) = readlines(abspath_(path))
loaddlm_(path::Nothing) = nothing
loaddlm_(path::AbstractString) = begin
	# empty strings is the same as providing nothing.
    path != "" || return nothing
    ReadWrite.loaddlm(abspath_(path); header=true)
end
loaddlm_(path::Nothing, T::Type) = nothing
loaddlm_(path::AbstractString, T::Type) = begin
	# empty strings is the same as providing nothing.
    path != "" || return nothing
    ReadWrite.loaddlm(abspath_(path), T; header=true)
end
loaddlm_(path::Nothing, colnames::Vector, rownames::Vector) = nothing
loaddlm_(path, colnames::Vector, rownames::Vector) = DataFrameUtils.subtable(loaddlm_(path), colnames, rownames)
savedlm_(path::AbstractString, x::AbstractArray; colnames=nothing, rownames=nothing) = begin
    ReadWrite.savedlm(abspath_(path), x; colnames=colnames, rownames=rownames)
end


"Get names of measured genes, genes mutated in each experiment and unique list of mutated genes."
function get_logFC_genes(logFC::DataFrame, delim)
    meas_genes = logFC[!, 1]
    @assert eltype(meas_genes) <: AbstractString "No row names provided in logFC file."
    
    exp_genes = names(logFC)[2:end]
    exp_genes = delim === nothing ? [[n] for n in exp_genes] : split.(exp_genes, delim)
    
    mut_genes = unique(n for ns in exp_genes for n in ns)

    meas_genes, exp_genes, mut_genes
end

"""
Match the column to the row names of X to make a logical matrix indicating 
- meas_genes: List of genes measured. Row names of X. 
- exp_genes: List of gene name sets that are mutated in each experiment. Corresponds to column names of X.
"""
function names2J(meas_genes::Vector{<:AbstractString}, exp_genes::Vector{Vector{T}} where T<:AbstractString)
    J = zeros(Bool, length(meas_genes), length(exp_genes))
    for i_exp in 1:length(exp_genes)
        for mutName in exp_genes[i_exp]
            J[meas_genes .== mutName, i_exp] .= true
        end
    end
    # There might be a spelling problem if there is no match for an experiment
    @assert all(sum(J, dims=1) .>= 1)
    J
end


"""
Read gene list, or match regex against a full list of genes.
"""
read_geneList(fnameORregex::Nothing, fullList::Vector{<:AbstractString}) = fullList
read_geneList(fnameORregex::String, fullList::Vector{<:AbstractString}) = begin
    isfile(fnameORregex) || return read_geneList(Regex(fnameORregex))
    readlines_(fnameORregex)
end
read_geneList(rx::Regex, fullList::Vector{<:AbstractString}) = begin
    matches = fullList[match.(rx, fullList) .!== nothing]
    @assert length(matches) > 0 "Nothing found for $rx."
    @info "$(length(matches)) genes matched on pattern $rx."
    return matches
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

"For tables with row names in first column."
table2mat(df::DataFrame, T::Type=Float64) = Matrix{T}(df[!, 2:end])


function infer(logFC::String, TF=nothing, KP=nothing, out_WT::String="WT_infer.mat", out_WP::String="WP_infer.mat";
        log::Union{<:IO,<:AbstractString}=stdout, mut_sep=nothing, epochs::Integer=5000, opt="ADAMW", 
        lr::Float64=0.001, decay::Real=0, WT=nothing, WP=nothing, WT_mask=nothing, WP_mask=nothing,
        lambda_Bstar::Real=.1, lambda_absW::Real=0., reg_WT::Bool=true, train_WT::Bool=true)
    
	logFC = loaddlm_(logFC)
    meas_genes, exp_genes, mut_genes = get_logFC_genes(logFC, mut_sep)
    TFs = read_geneList(TF, mut_genes)
    KPs = read_geneList(KP, mut_genes)
    @assert length(intersect(TFs, KPs)) == 0 "Genes as both TF and KP not implemented."
    TFKPs = [TFs; KPs]
    Os  = setdiff(meas_genes, TFKPs)
    TFKPOs = [TFKPs; Os]
    # W is easier to deal with if the diagonal is self edges so we sort TF, KP, O
    logFC = DataFrameUtils.rowselect(logFC, TFKPOs)
    
    X = table2mat(logFC)
    J = names2J(TFKPOs, exp_genes)
    @info "$(size(J, 1)) measured genes in $(size(J, 2)) experiments with a total of $(sum(J)) mutated genes."

    Wₜinit = loaddlm_(WT, TFKPOs, TFs)
    Wₚinit = loaddlm_(WP, TFKPs, KPs)
    Mₜ, Sₜ = mat2MS(loaddlm_(WT_mask, TFKPOs, TFs))
    Mₚ, Sₚ = mat2MS(loaddlm_(WP_mask, TFKPs, KPs))
    
    mdl = Model.get_model(length(TFs), length(KPs), length(Os), J; Wₜ=Wₜinit, Wₚ=Wₚinit, Mₜ=Mₜ, Mₚ=Mₚ, Sₜ=Sₜ, Sₚ=Sₚ)
    train_WT || Model.untrainable_Wₜ()
    # test that model works
    @assert size(mdl(X)) == size(X)
	
    opt = parse_optimizer(opt, lr, decay)

    Wₜ, Wₚ = GradientDescent.train(mdl, X, log; epochs=epochs, opt=opt, λBstar=lambda_Bstar, λabsW=lambda_absW, reg_Wₜ=reg_WT)
	
	train_WT && savedlm_(out_WT, Wₜ; colnames=TFs, rownames=TFKPOs)
	savedlm_(out_WP, Wₚ; colnames=KPs, rownames=TFKPs)
end

# parse args if run on command line as opposed to being imported
if abspath(PROGRAM_FILE) == @__FILE__
    ArgParseUtils.main(argument_parser, infer)
end

