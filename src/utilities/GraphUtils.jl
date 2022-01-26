#!/usr/bin/env julia
module GraphUtils
using Statistics
using DataFrames
using SparseArrays
using Colors: hex  # override hex so we can take hex of color
Main.@use "utilities/XGMML"
Main.@src "utilities/ColorUtils"
using ..ColorUtils: divergent_lerp
Main.@use "utilities/MathUtils"
Main.@src "simulation/GeneRegulation"
using ..GeneRegulation: estimate_Wₜ
Main.@use "inference/Model"
using ..Model: nₜnₚnₒ


const default_hspace, default_vspace = 45, 80
const bg_color = "#FFFFFF"

"""
Get a matrix containing text in node file format.
- nᵥ: total number of proteins
- nₜ: number of TFs
- nₚ: number of KPs
rows are identifiers.
"""
function nodes(nᵥ::Integer, nₜ::Integer, nₚ::Integer)
	labels = [["TF$i" for i in 1:nₜ]; ["KP$i" for i in 1:nₚ]; ["O$i" for i in 1:(nᵥ-(nₜ+nₚ))]]
	types = [["TF" for _ in 1:nₜ]; ["KP" for _ in 1:nₚ]; ["O" for _ in 1:(nᵥ-(nₜ+nₚ))]]
	DataFrame(label=labels, type=types)
end

function edges(matrix::SparseMatrixCSC)
	i,j,v = findnz(matrix)
	DataFrame(source=j, target=i, weight=v)
end
edges(matrix::Matrix) = edges(sparse(matrix))

"""
- space: minimum space that will be given horizontally.
"""
function xgmml_x(nₜ::Integer, nₚ::Integer, nₒ::Integer, space=default_hspace)
	width = space * max(nₜ, nₚ, nₒ)
	[
		[(i-.5) * width / nₜ for i in 1:nₜ];
		[(i-.5) * width / nₚ for i in 1:nₚ];
		[(i-.5) * width / nₒ for i in 1:nₒ]
	]
end
"""
- space: space that will be given vertically.
"""
xgmml_y(nₜ::Integer, nₚ::Integer, nₒ::Integer, space=default_vspace) = [
    [ 1space for i in 1:nₜ];
    [ 0space for i in 1:nₚ]; # placing KP above TF
    [ 2space for i in 1:nₒ]
]

function xgmml_labels(nₜ::Integer, nₚ::Integer, nₒ::Integer)
	pad = max(nₜ, nₚ, nₒ) |> string |> length
	[
	    ["TF"*lpad(i,pad,"0") for i in 1:nₜ];
	    ["KP"*lpad(i,pad,"0") for i in 1:nₚ];
	    ["O"*lpad(i,pad,"0") for i in 1:nₒ]
	]
end

function xgmml_fills(nₜ::Integer, nₚ::Integer, nₒ::Integer)
	[
	    ["#86ac32" for _ in 1:nₜ];
	    ["#af6fb6" for _ in 1:nₚ];
	    ["#e6dd47" for _ in 1:nₒ]
	]
end
"Get a hex color for each value in X where the color is a divergent color with limits min, max."
xgmml_fills(X::Matrix, min=minimum(X), max=maximum(X)) = "#" .* hex.(divergent_lerp.(X, min, max))

function xgmml_shapes(nₜ::Integer, nₚ::Integer, nₒ::Integer)
	[
	    ["ELLIPSE" for _ in 1:nₜ];
	    ["DIAMOND" for _ in 1:nₚ];
	    ["RECTANGLE" for _ in 1:nₒ]
	]
end


function xgmml_nodes(nₜ::Integer, nₚ::Integer, nₒ::Integer; x=xgmml_x(nₜ, nₚ, nₒ), y=xgmml_y(nₜ, nₚ, nₒ), labels=xgmml_labels(nₜ, nₚ, nₒ), fills=xgmml_fills(nₜ, nₚ, nₒ), shapes=xgmml_shapes(nₜ, nₚ, nₒ), extra_atts...)
	[XGMML.Node(x[i], y[i]; label=labels[i], fill=fills[i], shape=shapes[i], (k=>v[i] for (k,v) in extra_atts)...) for i in 1:nₜ+nₚ+nₒ]
end
xgmml_nodes(Wₜ, Wₚ; kwargs...) = xgmml_nodes(nₜnₚnₒ(Wₜ, Wₚ)...; kwargs...)
# allow WₜWₚ tuple so we can use net as either a Network or as (Wₜ, Wₚ) in the function xgmml(net, ...)
xgmml_nodes(WₜWₚ::Tuple; kwargs...) = xgmml_nodes(WₜWₚ...; kwargs...)
xgmml_nodes(net; kwargs...) = xgmml_nodes(net.nₜ, net.nₚ, net.nₒ;
    max_transcription=net.max_transcription, max_translation=net.max_translation,
    λ_mRNA=net.λ_mRNA, λ_prot=net.λ_prot, λ₊=[net.λ₊; zeros(net.nₒ)], λ₋=[net.λ₋; zeros(net.nₒ)], # assuming O is at the end of node list
    α₀=[g.α[1] for g in net.genes], labels=net.names, kwargs...)

edge_color(weight::Real) = "#" * hex(divergent_lerp(weight, -1, 1))
edge_color_TF(weight::Real) = edge_color(weight::Real)
edge_color_KP(weight::Real) = edge_color(weight::Real)
edge_color_TF(weight::AbstractString) = weight == "-" ? "#4B562D" : "#A1CC2B" 
edge_color_KP(weight::AbstractString) = weight == "-" ? "#6D3775" : "#AF70B6"
edge_opacity(weight::Real) = divergent(weight, -.1, .1) |> abs |> to256
edge_opacity(weight::AbstractString) = 255;
arrow_TF(weight::Real) = weight >= 0 ? "DELTA" : "T"
function arrow_TF(weight::AbstractString)
	if weight == "+" return "DELTA"
	elseif weight == "-" return "T"
	else return "DELTA" end # default value
end
arrow_KP(weight::Real) = weight >= 0 ? "CIRCLE" : "SQUARE"
function arrow_KP(weight::AbstractString)
	if weight == "+" return "CIRCLE"
	elseif weight == "-" return "SQUARE"
	else return "DELTA" end # default value
end

function xgmml_edges(Wₜ::Matrix, Wₚ::Matrix)
	nₜ = size(Wₜ, 2)

	TF_edges = [
        XGMML.Edge(source, target, bg_color; arrow=arrow_TF(weight), color=edge_color_TF(weight), opacity=edge_opacity(weight), weight=weight)
        for (target, source, weight) in zip(findnz(sparse(Wₜ))...)
    ]
	KP_edges = [
	    # + nₜ since the order of nodes should be TF, KP, O
        XGMML.Edge(source + nₜ, target, bg_color; arrow=arrow_KP(weight), color=edge_color_KP(weight), opacity=edge_opacity(weight), weight=weight)
        for (target, source, weight) in zip(findnz(sparse(Wₚ))...)
    ]

	[KP_edges; TF_edges]
end
# allow WₜWₚ tuple so we can use net as either a Network or as (Wₜ, Wₚ) in the function xgmml(net, ...)
xgmml_edges(WₜWₚ::Tuple) = xgmml_edges(WₜWₚ...)
xgmml_edges(net) = xgmml_edges(estimate_Wₜ(net), net.Wₚ₊ - net.Wₚ₋)

"""
Bend the edges if the target is 2 or more places away from the source within the same row of nodes.
Bend to the left relative to direction of edge, in order to avoid visual overlap when two nodes both have an edge to one another.
"""
function xgmml_bend!(graph::XGMML.Graph, bend=.3; hspace=default_hspace, vspace=default_vspace)
	for e in graph.edges
		source, target = XGMML.get_position(graph, e.source), XGMML.get_position(graph, e.target)
		seps = abs.(target .- source) ./ [hspace, vspace]  # seps = n spaces separating source and target
		if seps[1] >= 2 && seps[2] == 0
			XGMML.set_anchor(e, source, target, -bend)
		end
	end
end

# in order to allow xgmml take either (Wₜ, Wₚ) or Network.
# Has to be explicit with `Model.` to extend.
Model.nₜnₚnₒ(net::Tuple) = nₜnₚnₒ(net...)
Model.nₜnₚnₒ(net) = net.nₜ, net.nₚ, net.nₒ

"""
- net: either (Wₜ, Wₚ) or simulation Network
"""
function _graph(net; title="net")
	graph = XGMML.Graph(title, xgmml_nodes(net), xgmml_edges(net))
	xgmml_bend!(graph)
	graph
end
_graph(Wₜ::Matrix, Wₚ::Matrix; title="net") = _graph((Wₜ, Wₚ); title=title)

"Get a TF, KP, O network defined by its Wₜ and Wₚ in .xgmml format which can be imported into Cytoscape."
xgmml(Wₜ::Matrix, Wₚ::Matrix; title="net") = XGMML.xgmml(_graph(Wₜ, Wₚ; title=title))
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Nothing; title="net") = xgmml(Wₜ, Wₚ; title=title)
xgmml(net; title="net") = XGMML.xgmml(_graph(net; title=title))

"""
- net: either (Wₜ, Wₚ) or simulation Network
- X: each column is node values to visualize.
- highlight: index of node to highlight for each column in X.
"""
function xgmml(net, X::Matrix, highlight::Union{Nothing,Vector{Int}}=nothing; title="net")
	nₜ, nₚ, nₒ = nₜnₚnₒ(net)
	K = size(X, 2)
	fills = xgmml_fills(X, -1, 1)
	
	graphs::Vector{XGMML.Graph} = []
	for k in 1:K
		Δy = sum([nₜ,nₚ,nₒ] .> 0) * default_vspace * k
		nodes = xgmml_nodes(net, y=xgmml_y(nₜ, nₚ, nₒ) .+ Δy, fills=fills[:, k], value=X[:, k])
		edges = xgmml_edges(net)
        highlight === nothing || (nodes[highlight[k]].stroke_width = 5.)
		graph = XGMML.Graph(title, nodes, edges)
		xgmml_bend!(graph)
		push!(graphs, graph)
	end
	XGMML.xgmml(graphs)
end
xgmml(Wₜ::Matrix, Wₚ::Matrix, X::Matrix, highlight=nothing; title="net") = xgmml((Wₜ,Wₚ), X, highlight; title=title)


end;
