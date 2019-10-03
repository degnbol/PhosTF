#!/usr/bin/env julia
module Cytoscape
using DataFrames
using SparseArrays
using Colors: hex  # override hex so we can take hex of color
include("XGMML.jl"); using .XGMML: Graph
include("ColorUtils.jl"); using .ColorUtils: divergent_lerp
include("MathUtils.jl"); using .MathUtils

"""
Get a matrix containing text in node file format.
- n: total number of proteins
- nₜ: number of TFs
- nₚ: number of PKs
rows are identifiers.
"""
function nodes(n::Int, nₜ::Int, nₚ::Int)
	labels = [["TF$i" for i in 1:nₜ]; ["PK$i" for i in 1:nₚ]; ["X$i" for i in 1:(n-(nₜ+nₚ))]]
	types = [["TF" for _ in 1:nₜ]; ["PK" for _ in 1:nₚ]; ["X" for _ in 1:(n-(nₜ+nₚ))]]
	DataFrame(label=labels, type=types)
end

function edges(matrix::SparseMatrixCSC)
	i,j,v = findnz(matrix)
	DataFrame(source=j, target=i, weight=v)
end
edges(matrix::Matrix) = edges(sparse(matrix))


function xgmml_edges(Wₜ::Matrix, Wₚ::Matrix)
	nₜ = size(Wₜ,2)
	min_weight = minimum([minimum(Wₜ), minimum(Wₚ)])
	max_weight = maximum([maximum(Wₜ), maximum(Wₚ)])
	color(weight) = "#" * hex(divergent_lerp(weight, min_weight, max_weight))
	opacity(weight) = divergent(weight, min_weight, max_weight) |> abs |> to256
	TF_arrow(weight) = weight >= 0 ? "DELTA" : "T"
	PK_arrow(weight) = weight >= 0 ? "CIRCLE" : "SQUARE"
	
	id = 0
	TF_edges = [
	XGMML.Edge(source, target, id+=1, arrow=TF_arrow(weight), color=color(weight), opacity=opacity(weight))
	for (target,source,weight) in zip(findnz(sparse(Wₜ))...)]
	PK_edges = [
	XGMML.Edge(source+nₜ, target, id+=1, arrow=PK_arrow(weight), color=color(weight), opacity=opacity(weight))
	for (target,source,weight) in zip(findnz(sparse(Wₚ))...)]
	
	[TF_edges; PK_edges]
end


"Get a PK, TF, X network defined by its Wₜ and Wₚ in .xgmml format which can be imported into Cytoscape."
function xgmml(Wₜ::Matrix, Wₚ::Matrix)
	(n,nₜ), nₚ = size(Wₜ), size(Wₚ,2)
	pad = max(nₜ, nₚ, n-(nₜ+nₚ)) |> string |> length
	space = 100
	
	PKs = [XGMML.Node(i*space, 0space, nₜ+i, "PK"*lpad(i,pad,"0"), fill="#af6fb6", shape="DIAMOND") for i in 1:nₚ]
	TFs = [XGMML.Node(i*space, 1space, i, "TF"*lpad(i,pad,"0"), fill="#86ac32") for i in 1:nₜ]
	Xs  = [XGMML.Node(i*space, 2space, nₜ+nₚ+i, "X"*lpad(i,pad,"0"), fill="#e6dd47") for i in 1:n-(nₜ+nₚ)]
	
	graph = XGMML.Graph("pktfx", [TFs; PKs; Xs], xgmml_edges(Wₜ, Wₚ))
	XGMML.xgmml(graph)
end
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Nothing) = xgmml(Wₜ::Matrix, Wₚ::Matrix)

"""
- X: each column is node values to visualize (not to be confused with Xs referring to genes).
- highlight: a node to highlight for each column in X.
Set to nothing to not highlight anything. Default is to highlight the i-th protein in the i-th network.
The node is highlighted with a black border.
"""
function xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Matrix, highlight=1:size(X,2))
	(n,nₜ), nₚ, K = size(Wₜ), size(Wₚ,2), size(X,2)
	pad = max(nₜ, nₚ, n-(nₜ+nₚ)) |> string |> length
	space = 100
	pos(i,j) = (i*space, j*space)
	
	min_X = minimum(X)
	max_X = maximum(X)
	color(i,k) = "#" * hex(divergent_lerp(X[i,k], min_X, max_X))
	
	graphs::Vector{Graph} = []
	for k in 1:K
		PKs = [XGMML.Node(pos(i,0+3k)..., "PK"*lpad(i,pad,"0"), fill=color(nₜ+i,k), shape="DIAMOND") for i in 1:nₚ]
		TFs = [XGMML.Node(pos(i,1+3k)..., "TF"*lpad(i,pad,"0"), fill=color(i,k)) for i in 1:nₜ]
		Xs  = [XGMML.Node(pos(i,2+3k)...,  "X"*lpad(i,pad,"0"), fill=color(nₜ+nₚ+i,k)) for i in 1:n-(nₜ+nₚ)]
		graph = XGMML.Graph("pktfx", [TFs; PKs; Xs], xgmml_edges(Wₜ, Wₚ))
	end
	
	XGMML.xgmml(graphs)
end

end;
