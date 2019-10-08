#!/usr/bin/env julia
module FluxUtils
using LinearAlgebra
using Flux
import Base: *

export zerodiag
export random_weight

"""
A differentiable way to set the diagonal of a matrix to zeros.
"""
function zerodiag(matrix)
	remove = zeros(size(matrix))
	# collect converts the tracked values to normal values, meaning that we cannot differentiate this step.
	remove[diagind(remove)] = collect(diag(matrix))
	matrix - remove
end

"""
Set Gaussian random values where μ is set by the number of values.
"""
random_weight(n::Integer, m::Integer) = randn(n, m) ./ √(n*m)

# overloading missing multiplication signature that gives ambiguity.
# we are selecting the tracker multiplication over the LinearAlgebra one since we obviously want tracking.
*(x::LinearAlgebra.Diagonal, y::Tracker.TrackedArray{T,2,A} where A where T) = Tracker.track(*, x, y)

end;