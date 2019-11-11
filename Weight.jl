#!/usr/bin/env julia
if !isdefined(Main, :Model) include("Model.jl") end

"Construction of valid and random weight matrices, detection of silent edges and weight matrix correction."
module Weight
using LinearAlgebra
using SparseArrays
using ..Model

function random_Wₜ(n::Integer, nₜ::Integer, nₚ::Integer)
	mat = rand([-1, 0, 1], (n, nₜ))
	mat[diagind(mat, -nₚ)] .= 0
	mat
end

function random_Wₚ(nₜ::Integer, nₚ::Integer)
	mat = rand(Uniform(.5, 1.), (nₜ+nₚ, nₚ))
	mat[rand([false, true], size(mat))] .= 0
	mat[diagind(mat)] .= 0
	mat
end


self_loops(mat::AbstractMatrix, k::Integer=0) = (diag(mat, k) .!= 0)


"""
Mask of edges that leads to deadends, which are phosphorylation regulation leading to nodes that regulate nothing.
"""
function silent_edges(W::AbstractMatrix, nₚ::Integer)
	n = size(W,1)
	V = vec(all(W .== 0, dims=1))  # nodes not regulating anything
	E = falses(n,n)
	E_last = trues(n,n)  # just any value that is different from E
	while any(E .!= E_last)
		E_last .= E
		E[V,1:nₚ] .|= W[V,1:nₚ] .!= 0   # mark P edges as silent if they are onto a silent node
		V = vec(all((W .== 0) .| E, dims=1))  # recalculate silent nodes where the edges marked silent are removed from W
	end
	E
end

"""
Proteins regulated only by phosphatases, which means the phosphate edges serve no purpose since the protein will always be unphosphorylated.
return: 1D bit vector, length size(Wₚ,1). true for each protein regulated by only phosphatases (and at least 1), false otherwise.
"""
function silent_phosphatases(Wₚ::AbstractMatrix)
	kinase_reg, phosphate_reg = vec(any(Wₚ.>0, dims=2)), vec(any(Wₚ.<0, dims=2))
	phosphate_reg .& .!kinase_reg
end


"""
Check a calculation for issues and print a warning message if there are issues.
- arr: bool array, result from calculation of some test where true refers to an issue.
- msg: message to print in case of issues.
"""
function check(arr, msg)
	if any(arr)
		arr = findall(arr)
		@warn("$msg: $arr")
	end
end

"""
Set self-loops to zero and show a warning if self-loops are found.
return: bool indicating if any were found.
"""
function correct_self_loops!(W::AbstractMatrix)
	self_loops = self_loops(W)
	check(self_loops, "self loops in W")
	W[diagind(Wₚ)[self_loops]] .= 0
	return any(self_loops)
end
function correct_self_loops!(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	nₚ = size(Wₚ,2)
	self_loops_Wₜ, self_loops_Wₚ = self_loops(Wₜ,-nₚ), self_loops(Wₚ)
	check(self_loops_Wₜ, "self loops in Wₜ")
	check(self_loops_Wₚ, "self loops in Wₚ")
	Wₜ[diagind(Wₜ,-nₚ)[self_loops_Wₜ]] .= 0
	Wₚ[diagind(Wₚ)[self_loops_Wₚ]] .= 0
	return any(self_loops_Wₜ) || any(self_loops_Wₚ)
end
"""
Set silent (PK/PP) edges to zero and show a warning if any are found.
return: bool indicating if any were found.
"""
function correct_silent_edges!(W::AbstractMatrix, nₚ::Integer)
	E = silent_edges(W, nₚ)
	check(E, "silent edges")
	W[E] .= 0
	return any(E)
end
function correct_silent_edges!(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	_, nₜ, nₚ = Model.nₓnₜnₚ(Wₜ, Wₚ)
	W = Model._W(Wₜ, Wₚ)
	out = correct_silent_edges!(W, nₚ)
	Wₜ[:], Wₚ[:] = Model.WₜWₚ(W, nₜ, nₚ)
	return out
end
"""
Set (PK/PP) phosphate edges to zero if they will be silent in simulation and show a warning if any are found.
return: bool indicating if any were found.
"""
function correct_silent_phosphates!(Wₚ::AbstractMatrix)
	phosphates = silent_phosphatases(Wₚ)
	check(phosphates, "silent phosphate regulations")
	Wₚ[phosphates,:] .= 0
	return any(phosphates)
end
correct_silent_phosphates!(W::AbstractMatrix, nₚ::Integer) = correct_silent_phosphates!(View(W,:,1:nₚ))

"""
Correct self loops, silent edges, and silent phosphorylation regulations.
return: bool indicating if anything was corrected.
"""
function correct!(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	corrected = correct_self_loops!(Wₜ, Wₚ)
	while correct_silent_edges!(Wₜ, Wₚ) | correct_silent_phosphates!(Wₚ)
		corrected = true
	end
	corrected
end
function correct!(W::AbstractMatrix, nₚ::Integer)
	corrected = correct_self_loops!(W)
	while correct_silent_edges!(W, nₚ) | correct_silent_phosphates!(W, nₚ)
		corrected = true
	end
	corrected
end

threshold!(W::AbstractMatrix, thres::AbstractFloat=0.001) = W[abs.(W) .< thres] .= 0

end;
