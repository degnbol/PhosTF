#!/usr/bin/env julia
@src "simulation/GeneRegulation" # To get Network
@src "inference/Model"
@src "utilities/ReadWrite"
@src "utilities/ArrayUtils"
@src "utilities/GraphUtils"
using Fire
using DataFrames

"""
Write a graph defined by weight matrices to xgmml format.
See lower level function xgmml for arguments.
- Wₜ, Wₚ: fnames of tab-separated matrices with column and row names to indicate source and target, respectively.
"""
@main function xgmml(Wₜ::String, Wₚ::String; o=stdout, title=nothing, X=nothing, highlight=nothing, typecol=true, showlab=true)
	Wₜ, Wₚ = ReadWrite.loaddlm(Wₜ; header=true), ReadWrite.loaddlm(Wₚ; header=true)
    gene_names = Wₜ[!, 1]
    # minus one due to rownames
    nₜ, nₚ = ncol(Wₜ)-1, ncol(Wₚ)-1
    # check name order correspondence
    @assert all(gene_names[1:nₜ+nₚ] .== Wₚ[!, 1] .== [names(Wₜ)[2:end]; names(Wₚ)[2:end]])
    Wₜ = Matrix(Wₜ[:, 2:end])
    Wₚ = Matrix(Wₚ[:, 2:end])
	xgmml((Wₜ, Wₚ), gene_names, o, nₜ, nₚ, title, X, highlight, typecol ? nothing : "gray", showlab)
end
"""
- i: e.g. net.bson
See lower level function xgmml for arguments.
"""
@main function xgmml(i::String; o=stdout, title=nothing, X=nothing, highlight=nothing, typecol=true, showlab=true)
	net = ReadWrite.load(i, GeneRegulation.Network)
	xgmml(net, net.names, o, net.nₜ, net.nₚ, title, X, highlight, typecol ? nothing : "gray", showlab)
end

"""
Called by the two @main top-level functions.
- net: either (Wₜ, Wₚ) or Network.
- X: fname for tab-separated file with node values, e.g. log₂(fold-change) or steady state values.
    Each column of X is used for coloring a separate copy of the graph.
    Optionally the first column can contain row names, which should be node names. The rownames column name can be e.g. "_", "row", or "rownames"
- highlight: optional String name of node to highlight.
- fills: node fill color if X is not given. Default is coloring by type.
- showlab: show node labels (gene_names).
"""
function xgmml(net, gene_names, o, nₜ, nₚ, title=nothing, X=nothing, highlight=nothing, fills=nothing, showlab=true)
    if title === nothing
        title = o == stdout ? "PhosTF" : splitext(basename(o))[1]
    end
	if X === nothing
        write(o, GraphUtils.xgmml(net; title=title, fills=fills, labels=showlab ? nothing : []))
	else
		X = ReadWrite.loaddlm(X, Float64; header=true)
        if eltype(X[!, 1]) <: AbstractString
            rownames = X[!, 1]
            colnames = names(X)[2:end]
            typeof(net) == Tuple || @assert(all(rownames .== gene_names))
            # highlight a node for each column of X if they are referred to in the column names
            if highlight === nothing && colnames ⊆ rownames
                highlight = colnames
            end
            X = Matrix(X[:, 2:end])
        else
            # highlight each of the proteins if there are as many experiments as TFs+KPs
            if highlight === nothing && nₜ+nₚ == K
                highlight = (1:K)
            end
            X = Matrix(X)
        end
        highlight = get_highlight(highlight, size(X, 2), rownames)
        write(o, GraphUtils.xgmml(net, X, highlight; title=title, labels=showlab ? gene_names : []))
	end
end

# get the node indexes to highlight given various different allowed input forms (typeof(x)).
get_highlight(x::Nothing, K::Integer, names) = nothing
get_highlight(x::Integer, K::Integer, names) = [x for _ ∈ 1:K]
get_highlight(x::String, K::Integer, names::Vector{<:AbstractString}) = [only(findall(names .== x)) for _ ∈ 1:K]
get_highlight(x::String, K::Integer, names::Nothing) = error("Highlighting $x only possible if rownames are in X.")
get_highlight(x::Vector{String}, K::Integer, names::Vector{<:AbstractString}) = [only(findall(names .== n)) for n ∈ x]
get_highlight(x::Vector{String}, K::Integer, names::Nothing) = error("Highlighting $x only possible if rownames are in X.")
get_highlight(x::Union{UnitRange{Int},Vector{Int}}, K::Integer, names) = x


