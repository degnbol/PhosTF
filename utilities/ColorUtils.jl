#!/usr/bin/env julia
module ColorUtils
using Colors
# PerceptualColourMaps is using deprecated functions
using PerceptualColourMaps
include("MathUtils.jl"); using .MathUtils

export divergent_lerp


D7 = convert.(RGB, cmap("D7"))


function divergent_lerp(v, min=-1, max=1, center=0; cmap=D7)
	l = length(cmap)
	# set v to 1 when at min, half of cmap length when at center and 'end' when at max
	v = ((l - 1) * divergent(v, min, max, center) + (l + 1)) / 2
	cmap[clamp(Int(round(v)), 1, l)]
end


end;
