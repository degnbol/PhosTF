#!/usr/bin/env julia
module ColorUtils
using Colors
include("MathUtils.jl"); using .MathUtils
using DelimitedFiles: readdlm

export divergent_lerp
export lerp

# Build the divergent colormap named D7 with the package PerceptualColourMaps which at time of writing is using deprecated functions that cause warnings. 
# The colormap was compiled and stored in the file colormap_D7.txt from where it is faster to read it.
#= using PerceptualColourMaps # PerceptualColourMaps is using deprecated functions =#
#= D7 = convert.(RGB, cmap("D7")) =#
D7 = [RGB(rgb...) for rgb in eachcol(readdlm(joinpath(@__DIR__, "colormap_D7.txt"), '\t'))]

function divergent_lerp(v, min=-1, max=1, center=0; cmap=D7)
	l = length(cmap)
	# set v to 1 when at min, half of cmap length when at center and 'end' when at max
	v = ((l - 1) * divergent(v, min, max, center) + (l + 1)) / 2
	cmap[clamp(Int(round(v)), 1, l)]
end

"""
Mix two hex colors and return as hex.
- amount: 0 → hex1, 1 → hex2
"""
function lerp(hex1, hex2, amount::AbstractFloat)
	col = weighted_color_mean(amount, parse(Colorant, hex2), parse(Colorant, hex1))
	"#" * hex(col)
end
lerp(hex1, hex2, amount::Integer) = lerp(hex1, hex2, amount/255)

end;
