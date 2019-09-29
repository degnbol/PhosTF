#!/usr/bin/env julia
module ArrayUtils
using LinearAlgebra: diagind
using Random: shuffle
using Distributions

export binary, TruncNormal, eye

"""
Convert a decimal number to binary and return it as a 1D bool array.
Padded to achieve a given size.
"""
function binary(num, size)
	collect(string(num, base=2, pad=size)) .== '1'
end

"""
Adding some default values to μ and σ.
"""
TruncNormal(l, u) = TruncatedNormal((l+u)/2, (u-l)/6, l, u)

"""
Ability to create non-square identity matrices by using dense matrix type.
"""
function eye(n::Int, m::Int)
	mat = zeros(n, m)
	mat[diagind(mat)] .= 1
	mat
end
function eye(size::Tuple)
	mat = zeros(size)
	mat[diagind(mat)] .= 1
	mat
end
eye(matrix) = eye(size(matrix))

shuffle_rows(matrix) = matrix[shuffle(1:end), :]
shuffle_rows_(matrix) = @view matrix[shuffle(1:end), :]
shuffle_columns(matrix) = matrix[:, shuffle(1:end)]
shuffle_columns_(matrix) = @view matrix[:, shuffle(1:end)]


end;