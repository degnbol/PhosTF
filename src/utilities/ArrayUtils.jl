#!/usr/bin/env julia
module ArrayUtils
using LinearAlgebra: diagind
using Random: shuffle
using Distributions

export tological
export binary, TruncNormal, eye
export tostring
export hcatpad
export reorder

"""
Convert a numerical index to a logical index.
Examples:
tological([2, 1], 3) → [1, 1, 0]
tological([2], 3)    → [0, 1, 0]
tological( 2 , 3)    → [0, 1, 0]
"""
tological(indices, length) = (out = falses(length); out[indices] .= true; out)
tological(indices::Integer, length) = (out = falses(length); out[indices] = true; out)

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
function eye(n::Integer, m::Integer)
	mat = zeros(n, m)
	mat[diagind(mat)] .= 1
	mat
end
function eye(size::Tuple)
	mat = zeros(size)
	mat[diagind(mat)] .= 1
	mat
end
eye(matrix::AbstractMatrix) = eye(size(matrix))

shuffle_rows(matrix::AbstractMatrix) = matrix[shuffle(1:end), :]
shuffle_rows_(matrix::AbstractMatrix) = @view matrix[shuffle(1:end), :]
shuffle_columns(matrix::AbstractMatrix) = matrix[:, shuffle(1:end)]
shuffle_columns_(matrix::AbstractMatrix) = @view matrix[:, shuffle(1:end)]

"""
A conversion of array to string that writes empty array as [] instead of e.g. Any[].
"""
function tostring(array::AbstractArray)
	if isempty(array) return "[]" end
	string(array)
end

"""
hcat but instead of error when dimensions are mismatched a value is padded below short arrays.
"""
hcatpad(X::Base.Generator; pad=0) = hcatpad(collect(X); pad=pad)
function hcatpad(X::Vector; pad=0)
	n = maximum(size(x,1) for x in X)
	hcat([[x; fill(pad, n-size(x,1),size(x,2))] for x in X]...)
end
hcatpad(V::Vector{<:Number}; pad=0) = hcat(V)


"""
Reorder rows and columns in an adjacency matrix, making sure edge values belong to the same target and source.
Works for non-square adjacency matrix by ignoring indexes for rows or columns not in range. 
If some indexes are only found along the longest axis a warning is printed.
In other cases where some indexes were not found, an error is thrown.
Examples:
`reorder(rand(4,4), [3,4,1,2])` works by reordering rows and columns and all values are kept.
`reorder(rand(5,4), [3,4,1,2,5])` works and displays warning.
`reorder(rand(4,4), [3,4,1,2,5])` throws BoundsError.
`reorder(rand(4,4), [3,4,1])` works and removes a single row and column.
"""
function reorder(adj::AbstractMatrix, order)
	n, m = size(adj)
	row_idx = [i for i in order if i<=n]
	col_idx = [i for i in order if i<=m]
	if (length(row_idx) < length(order) && n>=m) || (length(col_idx) < length(order) && m>=n)
		throw(BoundsError("Indexing outside range of longest axis."))
	end
	if length(row_idx) < length(order) || length(col_idx) < length(order)
		@warn("Indexing outside range of shortest axis.")
	end
	adj[row_idx,:][:,col_idx]  # reorders rows then columns
end




end;