#!/usr/bin/env julia
"""
Module for code to export a network in .xgmml format.
"""
module XGMML
using Statistics: mean
import Base.show
import Base.merge
Main.@use "utilities/StringUtils"
Main.@use "utilities/ArrayUtils"
Main.@use "utilities/ColorUtils" # lerp

const default_h = 25.

mutable struct Node
	atts::Dict{Any,Any}
	graphics::Dict
	id
	label
	x
	y
	w
	h
	shape
	fill
	stroke_width
	stroke_color
	"- shape: ELLIPSE, DIAMOND, RECTANGLE, ..."
	function Node(x, y, attributes, graphics; id=nothing, label=id, fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=default_h, w=h, extra_atts...)
		new(merge(attributes, extra_atts), graphics, id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	end
	function Node(x, y, graphics; id=nothing, label=id, fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=default_h, w=h, extra_atts...)
		new(merge(attributes(id, label), extra_atts), graphics, id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	end
	function Node(x, y; id=nothing, label=id, fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=default_h, w=h, extra_atts...)
		new(merge(attributes(id, label), extra_atts), graphics(label), id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	end
	function Node(pos; id=nothing, label=id, fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=default_h, w=h, extra_atts...)
		new(merge(attributes(id, label), extra_atts), graphics(label), id, label, pos[1], pos[2], w, h, shape, fill, stroke, stroke_color)
	end
	
	attributes(shared_name=nothing, name=shared_name) = Dict("shared name" => shared_name, "name" => name)
	function graphics(label; label_color="#000000", select_color="#00FFFF", opacity=255, stroke_opacity=255, label_opacity=255, font="SansSerif,plain", font_size=12, stroke="SOLID", tooltip="")
		Dict(
		"NODE_LABEL" => label,
		"NODE_LABEL_FONT_FACE" => font,
		"NODE_LABEL_FONT_SIZE" => font_size,
		"NODE_LABEL_COLOR" => label_color,
		"NODE_SELECTED_PAINT" => select_color,
		"NODE_TRANSPARENCY" => opacity,
		"NODE_BORDER_TRANSPARENCY" => stroke_opacity,
		"NODE_LABEL_TRANSPARENCY" => label_opacity,
		"NODE_BORDER_STROKE" => stroke,
		"NODE_TOOLTIP" => tooltip
		)
	end
end


mutable struct Edge
	atts::Dict{Any,Any}
	graphics::Dict
	id
	label
	source
	target
	directed::Bool
	width
	color
	anchor
	function Edge(attributes, graphics, source, target; id=nothing, label="", directed=true, width=2., color="#000000", anchor=nothing, extra_atts...)
		new(merge(attributes, extra_atts), graphics, id, label, source, target, directed, width, color, anchor)
	end
	function Edge(graphics, source, target; id=nothing, label="", directed=true, width=2., color=nothing, anchor=nothing, extra_atts...)
		if color === nothing color = graphics["EDGE_TARGET_ARROW_UNSELECTED_PAINT"] end
		new(merge(attributes(id, label), extra_atts), graphics, id, label, source, target, directed, width, color, anchor)
	end
	function Edge(source, target; id=nothing, label="", directed=true, width=2., arrow="DELTA", color="#000000", opacity::Integer=255, anchor=nothing, extra_atts...)
		new(merge(attributes(id, label), extra_atts), graphics(label, arrow=arrow, color=color, opacity=opacity), id, label, source, target, directed, width, color, anchor)
	end
	"Provide background color in order to emulate transparency in arrow head colors so they match edge line color."
	function Edge(source, target, bg_color; id=nothing, label="", directed=true, width=2., arrow="DELTA", color="#000000", opacity::Integer=255, anchor=nothing, extra_atts...)
		new(merge(attributes(id, label), extra_atts), graphics(label, arrow=arrow, color=lerp(bg_color, color, opacity), opacity=opacity), id, label, source, target, directed, width, color, anchor)
	end
	
	function attributes(shared_name=nothing, name=shared_name; shared_interaction="", interaction=shared_interaction)
		Dict("shared name" => shared_name, "name" => name,
		"shared interaction" => shared_interaction, "interaction" => interaction)
	end
	
	"""
	- type: SOLID, LONG_DASH, ...
	- arrow: DELTA, T, CIRCLE, ...
	- color: unselected arrow head color
	"""
	function graphics(label="", tooltip=""; type="SOLID", arrow="DELTA", arrow_size=6., color="#000000", label_color="#000000", select_color="#00FFFF", opacity=255, label_opacity=255, font_size=10, font="Dialog,plain", curved=true, visible=true)
		Dict(
		"EDGE_LABEL" => label,
		"EDGE_TOOLTIP" => tooltip,
		"EDGE_LINE_TYPE" => type,
		"EDGE_TARGET_ARROW_SHAPE" => arrow,
		"EDGE_SOURCE_ARROW_SHAPE" => "NONE",
		"EDGE_TARGET_ARROW_SIZE" => arrow_size,
		"EDGE_SOURCE_ARROW_SIZE" => arrow_size,
		"EDGE_TARGET_ARROW_UNSELECTED_PAINT" => color,
		"EDGE_SOURCE_ARROW_UNSELECTED_PAINT" => color,
		"EDGE_STROKE_SELECTED_PAINT" => select_color,
		"EDGE_TARGET_ARROW_SELECTED_PAINT" => select_color,
		"EDGE_SOURCE_ARROW_SELECTED_PAINT" => select_color,
		"EDGE_LABEL_COLOR" => label_color,
		"EDGE_TRANSPARENCY" => opacity,
		"EDGE_LABEL_TRANSPARENCY" => label_opacity,
		"EDGE_LABEL_FONT_FACE" => font,
		"EDGE_LABEL_FONT_SIZE" => font_size,
		"EDGE_CURVED" => curved,
		"EDGE_VISIBLE" => visible
		)
	end
end


struct Graph
	nodes::Vector{Node}
	edges::Vector{Edge}
	atts::Dict
	graphics::Dict
	name::String
	directed::Bool
	function Graph(name::String, nodes::Vector{Node}, edges::Vector{Edge}, bg_color="#FFFFFF"; title=name)
		atts = Dict("shared name" => name, "name" => name, "__Annotations" => [])
		graphics = Dict("NETWORK_TITLE" => title, "NETWORK_BACKGROUND_PAINT" => bg_color)
		set_ids(nodes)
		set_ids(edges)
		new(nodes, edges, atts, graphics, name, true)
	end
end

begin # printing
	att_repr(name, value) = """<att name="$name" value="$value" type="string" cy:type="String"/>\n"""
	function att_repr(name, value::Array)
		value = ArrayUtils.tostring(value)
		"""<att name="$name" value="$value" type="list" cy:type="List" cy:elementType="String"/>\n"""
	end
	att_repr(atts::Dict) = prod(att_repr(name,value) for (name,value) in atts)
	function Base.show(io::IO, n::Node)
		graphics = """<graphics fill="$(n.fill)" x="$(n.x)" y="$(n.y)" w="$(n.w)" h="$(n.h)" type="$(n.shape)" width="$(n.stroke_width)" outline="$(n.stroke_color)">\n""" * indent(att_repr(n.graphics)) * "</graphics>\n"
		
		print(io, """<node id="$(n.id)" label="$(n.label)">\n""" * indent(att_repr(n.atts) * graphics) * "</node>\n")
	end
	function edge_graphics_repr(e::Edge)
		out = """<graphics width="$(e.width)" fill="$(e.color)">
		""" * indent(att_repr(e.graphics))
		if e.anchor !== nothing out *= indent(att_repr("EDGE_BEND", join(e.anchor,','))) end
		out *= "</graphics>\n"
		out
	end
	function Base.show(io::IO, e::Edge)
		directed = e.directed ? 1 : 0
		print(io, """<edge id="$(e.id)" label="$(e.label)" source="$(e.source)" target="$(e.target)" cy:directed="$directed">\n""" * indent(att_repr(e.atts) * edge_graphics_repr(e)) * "</edge>\n")
	end
	function Base.show(io::IO, g::Graph)
		directed = g.directed ? 1 : 0
		header = """<graph label="$(g.name)" directed="$directed" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">
		"""
		graphics = "<graphics>\n" * indent(att_repr(g.graphics)) * "</graphics>\n"
		nodes = prod(repr(node) for node in g.nodes)
		edges = prod(repr(edge) for edge in g.edges)
		
		print(io, header * indent(att_repr(g.atts) * graphics * nodes * edges) * "</graph>\n")
	end
end


begin # getters
	get_node_ids(g::Graph) = [node.id for node in g.nodes]
	get_node_ids(gs::Vector{Graph}) = [node.id for g in gs for node in g.nodes]
	get_edge_ids(g::Graph) = [edge.id for edge in g.edges]
	get_edge_ids(gs::Vector{Graph}) = [edge.id for g in gs for edge in g.edges]
	
	"Check if there are duplicate ids among graphs"
	function has_duplicate_ids(g::Union{Graph,Vector{Graph}})
		node_ids = get_node_ids(g)
		edge_ids = get_edge_ids(g)
		length(node_ids) > length(unique(node_ids)) || length(edge_ids) > length(unique(edge_ids))
	end
	
	"Center of mass."
	get_center(ns::Vector{Node}) = [mean(n.x for n in ns), mean(n.y for n in ns)]
	get_position(g::Graph, node::Integer) = [g.nodes[node].x, g.nodes[node].y]
end


begin # setters
	"If nodes/edges are without id, set it to their position in the vector."
	function set_ids(nes::Union{Vector{Node},Vector{Edge}})
		for i in 1:length(nes)
			if nes[i].id === nothing set_id(nes[i], i) end
		end
	end
	function set_id(ne, id)
		ne.id = id
		if ne.label === nothing set_label(ne, id) end
		if ne.atts["shared name"] === nothing set_shared_name(ne, id) end
	end
	function set_label(n::Node, label)
		n.label = label
		if n.graphics["NODE_LABEL"] === nothing n.graphics["NODE_LABEL"] = label end
	end
	function set_label(e::Edge, label)
		e.label = label
		if e.graphics["EDGE_LABEL"] === nothing e.graphics["EDGE_LABEL"] = label end
	end
	function set_shared_name(ne, name)
		ne.atts["shared name"] = name
		if ne.atts["name"] === nothing ne.atts["name"] = name end
	end
	set_shared_name(ne, name::Integer) = set_shared_name(ne, string(name))

	translate_x(g::Graph, Δ) = for node in g.nodes translate_x(node, Δ) end
	translate_y(g::Graph, Δ) = for node in g.nodes translate_y(node, Δ) end
	translate_x(n::Node, Δ) = n.x += Δ
	translate_y(n::Node, Δ) = n.y += Δ
	increment_ids(es::Vector{Edge}, Δ::Int) = for e in es e.id += Δ end
	append_ids(es::Vector{Edge}, Δ::String) = for e in es e.id = "$(e.id)$Δ" end
	function increment_node_ids(g::Graph, Δ::Int)
		for n in g.nodes n.id += Δ end
		# make sure edges refer to the correct source and target post increment.
		increment_source(g.edges, Δ)
		increment_target(g.edges, Δ)
	end
	function append_node_ids(g::Graph, Δ::String)
		for n in g.nodes n.id = "$(n.id)$Δ" end
		# make sure edges refer to the correct source and target post change.
		append_source(g.edges, Δ)
		append_target(g.edges, Δ)
	end
	increment_source(es::Vector{Edge}, Δ::Integer) = for e in es e.source += Δ end
	increment_target(es::Vector{Edge}, Δ::Integer) = for e in es e.target += Δ end
	append_source(es::Vector{Edge}, Δ::String) = for e in es e.source = "$(e.source)$Δ" end
	append_target(es::Vector{Edge}, Δ::String) = for e in es e.target = "$(e.target)$Δ" end
	
	function set_anchor(e::Edge, source::Vector, target::Vector, bend::Real)
		v = target - source
		v̂ = [-v[2], v[1]] # hat the vector. 90 degree rotation
		e.anchor = @. source + .5v + bend*v̂
	end
	set_anchor(e::Edge, source::Node, target::Node, bend::Real) = set_anchor(e, [source.x, source.y], [target.x, target.y], bend)
end


begin # methods
	"""
	Merge graphs to one.
	Makes sure ids will not conflict between graphs.
	Displays non-terminating error if graphs are not identical except in their attributes etc.
	Uses first graph to append others to.
	"""
	function merge(graphs::Vector{Graph})
		if !are_similar(graphs) @error("Graphs are not similar") end
		if !has_duplicate_ids(graphs) return _merge(graphs) end
		resolve(graphs)
		_merge(graphs)
	end
	"""
	Merge graphs by adding nodes and edges to the first graph in the list.
	Does NOT check id duplicates.
	"""
	function _merge(graphs::Vector{Graph})
		for g in graphs[2:end]
			append!(graphs[1].nodes, g.nodes)
			append!(graphs[1].edges, g.edges)
		end
		graphs[1]
	end
	"Return if the graphs are identical in attributes."
	function are_similar(g1::Graph, g2::Graph)
		[g1.atts, g1.graphics, g1.directed] == [g2.atts, g2.graphics, g2.directed]
	end
	function are_similar(graphs::Vector{Graph})
		for g in graphs[2:end]
			if !are_similar(graphs[1], g) return false end
		end
		true
	end
	"""
	Resolve id conflicts.
	If conflicting ids are integers, they will be incremented by the sizes of the preceeding graphs in the list.
	Otherwise "_<graph number>" is appended.
	Does not check x, y.
	"""
	resolve(graphs::Vector{Graph}) = resolve(graphs, get_node_ids(graphs), get_edge_ids(graphs))
	"Int resolution."
	function resolve(graphs::Vector{Graph}, node_ids::Vector{Int}, edge_ids::Vector{Int})
		node_increment = length(graphs[1].nodes)
		edge_increment = length(graphs[1].edges)
		for graph in graphs[2:end]
			increment_node_ids(graph, node_increment)
			increment_ids(graph.edges, edge_increment)
			node_increment += length(graph.nodes)
			edge_increment += length(graph.edges)
		end
	end
	"String resolution"
	function resolve(graphs::Vector{Graph}, node_ids, edge_ids)
		for i in 1:length(graphs)
			append_node_ids(graphs[i], "_$i")
			append_ids(graphs[i].edges, "_$i")
		end
	end

	xgmml(graph::Graph) = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n""" * repr(graph)
	xgmml(graphs::Vector{Graph}) = xgmml(merge(graphs))





end

end;
