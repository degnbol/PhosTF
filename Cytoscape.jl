#!/usr/bin/env julia
include("utilities/XGMML.jl")
include("utilities/ColorUtils.jl")
include("utilities/MathUtils.jl")
if !isdefined(Main, :GeneRegulation) include("simulation/GeneRegulation.jl") end

module Cytoscape
using Statistics
using DataFrames
using SparseArrays
using Colors: hex  # override hex so we can take hex of color
import ..XGMML
using ..ColorUtils: divergent_lerp
using ..MathUtils
using ..GeneRegulation: nₓnₜnₚ, estimate_Wₜ


const default_space = 100

"""
Get a matrix containing text in node file format.
- n: total number of proteins
- nₜ: number of TFs
- nₚ: number of PKs
rows are identifiers.
"""
function nodes(n::Integer, nₜ::Integer, nₚ::Integer)
	labels = [["TF$i" for i in 1:nₜ]; ["PK$i" for i in 1:nₚ]; ["X$i" for i in 1:(n-(nₜ+nₚ))]]
	types = [["TF" for _ in 1:nₜ]; ["PK" for _ in 1:nₚ]; ["X" for _ in 1:(n-(nₜ+nₚ))]]
	DataFrame(label=labels, type=types)
end

function edges(matrix::SparseMatrixCSC)
	i,j,v = findnz(matrix)
	DataFrame(source=j, target=i, weight=v)
end
edges(matrix::Matrix) = edges(sparse(matrix))


xgmml_x(nₓ::Integer, nₜ::Integer, nₚ::Integer, space=default_space) = [[i*space for i in 1:nₚ]; [i*space for i in 1:nₜ]; [i*space for i in 1:nₓ]]
xgmml_y(nₓ::Integer, nₜ::Integer, nₚ::Integer, space=default_space) = [[ 0space for i in 1:nₚ]; [ 1space for i in 1:nₜ]; [ 2space for i in 1:nₓ]]

function xgmml_labels(nₓ::Integer, nₜ::Integer, nₚ::Integer)
	pad = max(nₓ, nₜ, nₚ) |> string |> length
	[["PK"*lpad(i,pad,"0") for i in 1:nₚ];
	 ["TF"*lpad(i,pad,"0") for i in 1:nₜ];
	 ["X" *lpad(i,pad,"0") for i in 1:nₓ]]
end

function xgmml_fills(nₓ::Integer, nₜ::Integer, nₚ::Integer)
	[["#af6fb6" for _ in 1:nₚ];
	 ["#86ac32" for _ in 1:nₜ];
	 ["#e6dd47" for _ in 1:nₓ]]
end
"Get a hex color for each value in X where the color is a divergent color with limits min, max."
xgmml_fills(X::Matrix, min=minimum(X), max=maximum(X)) = "#" .* hex.(divergent_lerp.(X, min, max))

function xgmml_shapes(nₓ::Integer, nₜ::Integer, nₚ::Integer)
	[["DIAMOND" for _ in 1:nₚ];
	 ["ELLIPSE" for _ in 1:nₜ];
	 ["RECTANGLE" for _ in 1:nₓ]]
end


function xgmml_nodes(nₓ::Integer, nₜ::Integer, nₚ::Integer; x=xgmml_x(nₓ, nₜ, nₚ), y=xgmml_y(nₓ, nₜ, nₚ), labels=xgmml_labels(nₓ, nₜ, nₚ), fills=xgmml_fills(nₓ, nₜ, nₚ), shapes=xgmml_shapes(nₓ, nₜ, nₚ), extra_atts...)
	nodes::Vector{XGMML.Node} = []
	[XGMML.Node(x[i], y[i]; label=labels[i], fill=fills[i], shape=shapes[i], (k=>v[i] for (k,v) in extra_atts)...)
	for i in 1:nₓ+nₚ+nₜ]
end
xgmml_nodes(WₜWₚ::Array; kwargs...) = xgmml_nodes(nₓnₜnₚ(WₜWₚ)...; kwargs...)
xgmml_nodes(WₜWₚ::Tuple; kwargs...) = xgmml_nodes(nₓnₜnₚ(WₜWₚ)...; kwargs...)
function xgmml_nodes(net; kwargs...)
	xgmml_nodes(nₓnₜnₚ(net)...;
	max_transcription=net.max_transcription, max_translation=net.max_translation,
	λ_mRNA=net.λ_mRNA, λ_prot=net.λ_prot, λ_phos=[net.λ_phos; zeros(net.nₓ)],
	phos_activation=[net.phos_activation; fill(missing, net.nₓ)], α₀=[g.α[1] for g in net.genes], kwargs...)
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
xgmml_edges(WₜWₚ::Array) = xgmml_edges(WₜWₚ...)
xgmml_edges(WₜWₚ::Tuple) = xgmml_edges(WₜWₚ...)
xgmml_edges(net) = xgmml_edges(estimate_Wₜ(net), net.Wₚₖ-net.Wₚₚ)

"""
Bend the edges if the target is 2 or more places away from the source within the same row of nodes.
Bend to the left relative to direction of edge, in order to avoid visual overlap when two nodes both have an edge to one another.
"""
function xgmml_bend!(graph::XGMML.Graph, bend=.3; space=default_space)
	for e in graph.edges
		source, target = XGMML.get_position(graph, e.source), XGMML.get_position(graph, e.target)
		seps = abs.(target .- source) ./ space  # seps = n spaces separating source and target
		if seps[1] >= 2 && seps[2] == 0
			XGMML.set_anchor(e, source, target, -bend)
		end
	end
end

function _graph(net; title="net")
	graph = XGMML.Graph(title, xgmml_nodes(net), xgmml_edges(net))
	xgmml_bend!(graph)
	graph
end
_graph(Wₜ::Matrix, Wₚ::Matrix; title="net") = _graph([Wₜ, Wₚ]; title=title)

"Get a PK, TF, X network defined by its Wₜ and Wₚ in .xgmml format which can be imported into Cytoscape."
xgmml(Wₜ::Matrix, Wₚ::Matrix; title="net") = XGMML.xgmml(_graph(Wₜ, Wₚ; title=title))
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Nothing; title="net") = xgmml(Wₜ, Wₚ; title=title)
xgmml(net; title="net") = XGMML.xgmml(_graph(net; title=title))

"""
- net: [Wₜ,Wₚ] or simulation Network
- X: each column is node values to visualize (not to be confused with Xs referring to genes).
- highlight: index(es) of node(s) to highlight for each column in X.
"""
function xgmml(net, X::Matrix, highlight=nothing; title="net")
	nₓ, nₜ, nₚ = nₓnₜnₚ(net); K = size(X,2)
	fills = xgmml_fills(X, -1, 1)
	
	graphs::Vector{XGMML.Graph} = []
	for k in 1:K
		Δy = sum([nₓ,nₜ,nₚ] .> 0) * default_space*k
		nodes = xgmml_nodes(net, y=xgmml_y(nₓ,nₜ,nₚ) .+ Δy, fills=fills[:,k], value=X[:,k])
		edges = xgmml_edges(net)
		if highlight != nothing nodes[highlight[k]].stroke_width = 5. end
		graph = XGMML.Graph(title, nodes, edges)
		xgmml_bend!(graph)
		push!(graphs, graph)
	end
	XGMML.xgmml(graphs)
end
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Matrix, highlight=nothing; title="net") = xgmml((Wₜ,Wₚ), X, highlight; title=title)

end;
