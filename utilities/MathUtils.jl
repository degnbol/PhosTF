#!/usr/bin/env julia
module MathUtils

export divergent
export to256

"get -1 when v is at min, 0 when it is at center and 1 when at max"
divergent(v, min=-1, max=1, center=0) = (v - center) / (v >= center ? (max - center) : (center - min))

"""
Convert float ∈ [0,1] to int in ∈ [0,255]
"""
to256(v) = clamp(Int(round(255v)), 0, 255)

end;
