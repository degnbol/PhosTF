#!/usr/bin/env julia
"""
Module for code to export a network in .xgmml format.
"""
module XGMML
import Base.show
include("StringUtils.jl"); using .StringUtils
include("ArrayUtils.jl")

"Attribute that can be added to a network, node or edge."
struct Attributes
	atts::Dict
	Attributes(ps::Pair{K,V}...) where {K,V} = new(Dict(ps))
	Attributes(ps::Pair...) = new(Dict(ps))
	
	repr(name, value) = """<att name="$name" value="$value" type="string" cy:type="String"/>\n"""
	function repr(name, value::Array)
		value = ArrayUtils.tostring(value)
		"""<att name="$name" value="$value" type="list" cy:type="List" cy:elementType="String"/>\n"""
	end
	function Base.show(io::IO, a::Attributes)
		print(io, prod(repr(name, value) for (name,value) in a.atts))
	end
end


struct Node
	atts::Attributes
	graphics::Attributes
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
	Node(x, y, attributes, graphics, id, label=id; fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=35.0, w=h) = new(attributes, graphics, id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	Node(x, y, graphics, id, label=id; fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=35.0, w=h) = new(attributes(id), graphics, id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	Node(x, y, id, label=id; fill="#999999", shape="ELLIPSE", stroke=0., stroke_color="#000000", h=35.0, w=h) = new(attributes(id), graphics(label), id, label, x, y, w, h, shape, fill, stroke, stroke_color)
	
	attributes(shared_name, name=shared_name) = Attributes("shared_name" => shared_name, "name" => name)
	function graphics(label; label_color="#000000", select_color="#FFFF00", opacity=255, stroke_opacity=255, label_opacity=255, font="SansSerif,plain", font_size=12, stroke="SOLID", tooltip="")
		Attributes(
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
	
	function Base.show(io::IO, n::Node)
		graphics = """<graphics fill="$(n.fill)" x="$(n.x)" y="$(n.y)" w="$(n.w)" h="$(n.h)" type="$(n.shape)" width="$(n.stroke_width)" outline="$(n.stroke_color)">\n""" * indent(repr(n.graphics)) * "</graphics>\n"
		
		print(io, """<node id="$(n.id)" label="$(n.label)">\n""" * indent(repr(n.atts) * graphics) * "</node>\n")
	end
end


struct Edge
	atts::Attributes
	graphics::Attributes
	id
	label
	source
	target
	directed::Bool
	width
	color
	Edge(attributes, graphics, source, target, id; label="", directed=true, width=2., color="#000000") = new(attributes, graphics, id, label, source, target, directed, width, color)
	function Edge(graphics, source, target, id; label="", directed=true, width=2., color=nothing)
		if color == nothing color = graphics["EDGE_TARGET_ARROW_UNSELECTED_PAINT"] end
		new(attributes(id, label), graphics, id, label, source, target, directed, width, color)
	end
	Edge(source, target, id; label="", directed=true, width=2., arrow="DELTA", color="#000000", opacity=255) = new(attributes(id, label), graphics(label, arrow=arrow, color=color, opacity=opacity), id, label, source, target, directed, width, color)
	
	function attributes(shared_name, name=shared_name; shared_interaction="", interaction=shared_interaction)
		Attributes(
		"shared name" => shared_name,
		"name" => name,
		"shared interaction" => shared_interaction,
		"interaction" => interaction
		)
	end
	
	"""
	- type: SOLID, LONG_DASH, ...
	- arrow: DELTA, T, CIRCLE, ...
	"""
	function graphics(label="", tooltip=""; type="SOLID", arrow="DELTA", arrow_size=6., color="#000000", label_color="#000000", select_color="#FFFF00", opacity=255, label_opacity=255, font_size=10, font="Dialog,plain", curved=true, visible=true)
		Attributes(
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
	
	function Base.show(io::IO, e::Edge)
		graphics = """<graphics width="$(e.width)" fill="$(e.color)">
		""" * indent(repr(e.graphics)) * "</graphics>\n"
		directed = e.directed ? 1 : 0
		
		print(io, """<edge id="$(e.id)" label="$(e.label)" source="$(e.source)" target="$(e.target)" cy:directed="$directed">\n""" * indent(repr(e.atts) * graphics) * "</edge>\n")
	end
end


struct Graph
	nodes::Vector{Node}
	edges::Vector{Edge}
	atts::Attributes
	graphics::Attributes
	name::String
	directed::Bool
	function Graph(name::String, nodes::Vector{Node}, edges::Vector{Edge}, bg_color="#FFFFFF"; title=name)
		atts = Attributes("shared name" => name, "name" => name, "__Annotations" => [])
		graphics = Attributes("NETWORK_TITLE" => title, "NETWORK_BACKGROUND_PAINT" => bg_color)
		new(nodes, edges, atts, graphics, name, true)
	end
	
	function Base.show(io::IO, g::Graph)
		directed = g.directed ? 1 : 0
		header = """<graph label="$(g.name)" directed="$directed" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">
		"""
		graphics = "<graphics>\n" * indent(repr(g.graphics)) * "</graphics>\n"
		nodes = prod(repr(node) for node in g.nodes)
		edges = prod(repr(edge) for edge in g.edges)
		
		print(io, header * indent(repr(g.atts) * graphics * nodes * edges) * "</graph>\n")
	end
end


xgmml(graph::Graph) = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n""" * repr(graph)

end;
