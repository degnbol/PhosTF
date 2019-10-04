#!/usr/bin/env julia
module Cytoscape
using Statistics
using DataFrames
using SparseArrays
using Colors: hex  # override hex so we can take hex of color
include("utilities/XGMML.jl"); using .XGMML: Graph
include("utilities/ColorUtils.jl"); using .ColorUtils: divergent_lerp
include("utilities/MathUtils.jl"); using .MathUtils
include("Model.jl"); using .Model: nₓnₜnₚ

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


xgmml_x(nₓ::Int, nₜ::Int, nₚ::Int, space=100) = [[i*space for i in 1:nₚ]; [i*space for i in 1:nₜ]; [i*space for i in 1:nₓ]]
xgmml_y(nₓ::Int, nₜ::Int, nₚ::Int, space=100) = [[ 0space for i in 1:nₚ]; [ 1space for i in 1:nₜ]; [ 2space for i in 1:nₓ]]

function xgmml_labels(nₓ::Int, nₜ::Int, nₚ::Int)
	pad = max(nₓ, nₜ, nₚ) |> string |> length
	[["PK"*lpad(i,pad,"0") for i in 1:nₚ];
	 ["TF"*lpad(i,pad,"0") for i in 1:nₜ];
	 ["X" *lpad(i,pad,"0") for i in 1:nₓ]]
end

function xgmml_fills(nₓ::Int, nₜ::Int, nₚ::Int)
	[["#af6fb6" for _ in 1:nₚ];
	 ["#86ac32" for _ in 1:nₜ];
	 ["#e6dd47" for _ in 1:nₓ]]
end
"Get a hex color for each value in X where the color is a divergent color with limits min, max."
xgmml_fills(X::Matrix, min=minimum(X), max=maximum(X)) = "#" .* hex.(divergent_lerp.(X, min, max))

function xgmml_shapes(nₓ::Int, nₜ::Int, nₚ::Int)
	[["DIAMOND" for _ in 1:nₚ];
	 ["ELLIPSE" for _ in 1:nₜ];
	 ["RECTANGLE" for _ in 1:nₓ]]
end


function xgmml_nodes(nₓ::Int, nₜ::Int, nₚ::Int; x=xgmml_x(nₓ, nₜ, nₚ), y=xgmml_y(nₓ, nₜ, nₚ), labels=xgmml_labels(nₓ, nₜ, nₚ), fills=xgmml_fills(nₓ, nₜ, nₚ), shapes=xgmml_shapes(nₓ, nₜ, nₚ), extra_atts...)
	nodes::Vector{XGMML.Node} = []
	[XGMML.Node(x[i], y[i]; label=labels[i], fill=fills[i], shape=shapes[i], (k=>v[i] for (k,v) in extra_atts)...)
	for i in 1:nₓ+nₚ+nₜ]
end

function xgmml_edges(Wₜ::Matrix, Wₚ::Matrix)
	nₚ = size(Wₚ,2)
	min_weight = minimum([minimum(Wₜ), minimum(Wₚ)])
	max_weight = maximum([maximum(Wₜ), maximum(Wₚ)])
	color(weight) = "#" * hex(divergent_lerp(weight, min_weight, max_weight))
	opacity(weight) = divergent(weight, min_weight, max_weight) |> abs |> to256
	TF_arrow(weight) = weight >= 0 ? "DELTA" : "T"
	PK_arrow(weight) = weight >= 0 ? "CIRCLE" : "SQUARE"
	
	PK_edges = [
	XGMML.Edge(source, target, arrow=PK_arrow(weight), color=color(weight), opacity=opacity(weight), weight=weight)
	for (target,source,weight) in zip(findnz(sparse(Wₚ))...)]
	TF_edges = [
	XGMML.Edge(source+nₚ, target, arrow=TF_arrow(weight), color=color(weight), opacity=opacity(weight), weight=weight)
	for (target,source,weight) in zip(findnz(sparse(Wₜ))...)]
	
	[PK_edges; TF_edges]
end


"Get a PK, TF, X network defined by its Wₜ and Wₚ in .xgmml format which can be imported into Cytoscape."
function xgmml(Wₜ::Matrix, Wₚ::Matrix)
	XGMML.xgmml(XGMML.Graph("pktfx", xgmml_nodes(nₓnₜnₚ(Wₜ, Wₚ)...), xgmml_edges(Wₜ, Wₚ)))
end
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Nothing) = xgmml(Wₜ::Matrix, Wₚ::Matrix)

"""
- X: each column is node values to visualize (not to be confused with Xs referring to genes).
- highlight: index(es) of node(s) to highlight for each column in X.
"""
function xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Matrix, highlight=nothing)
	nₓ, nₜ, nₚ = nₓnₜnₚ(Wₜ, Wₚ); K = size(X,2)
	fills = xgmml_fills(X, -1, 1)
	
	graphs::Vector{Graph} = []
	for k in 1:K
		Δy = sum([nₓ,nₜ,nₚ] .> 0) * 100k
		nodes = xgmml_nodes(nₓ,nₜ,nₚ, y=xgmml_y(nₜ,nₚ,nₓ) .+ Δy, fills=fills[:,k], value=X[:,k])
		if highlight != nothing nodes[highlight[k]].stroke_width = 5. end
		graph = XGMML.Graph("pktfx", nodes, xgmml_edges(Wₜ, Wₚ))
		push!(graphs, graph)
	end
	XGMML.xgmml(graphs)
end

end;
