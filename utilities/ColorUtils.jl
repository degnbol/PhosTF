#!/usr/bin/env julia
module ColorUtils
using Colors
using Suppressor
# PerceptualColourMaps is using deprecated functions
@suppress using PerceptualColourMaps
include("MathUtils.jl"); using .MathUtils

export divergent_lerp


D7 = convert.(RGB, cmap("D7"))


function divergent_lerp(v, min=-1, max=1, center=0; cmap=D7)
	l = length(cmap)
	# set v to 1 when -1, half of cmap length when 0 and 'end' when 1
	v = ((l - 1) * divergent(v) + (l + 1)) / 2
	cmap[clamp(Int(round(v)), 1, l)]
end


end;
