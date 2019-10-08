#!/usr/bin/env julia
"Construction of valid and random weight matrices, detection of silent edges and weight matrix correction."
module Weight
using LinearAlgebra
using SparseArrays
include("Model.jl"); using .Model: nₓnₜnₚ

function Wₜ(n::Integer, nₜ::Integer, nₚ::Integer)
	mat = rand([-1, 0, 1], (n, nₜ))
	mat[diagind(mat, -nₚ)] .= 0
	mat
end

function Wₚ(nₜ::Integer, nₚ::Integer)
	mat = rand(Uniform(.5, 1.), (nₜ+nₚ, nₚ))
	mat[rand([false, true], size(mat))] .= 0
	mat[diagind(mat)] .= 0
	mat
end


self_loops(mat, k::Integer=0) = (diag(mat, k) .!= 0)

"""
Return silent edges.
An edge is silent if it is part of a cascade not leading to any gene expression, i.e. it makes no difference on gene expression if the edge is present or not.
TF edges are never silent.
- Iₜ: diagonal matrix with 1s for each TF and 0s for any other proteins.
- Iₚ: diagonal matrix with 1s for each PK/PP and 0s for any other proteins.
"""
function silent_edges(W, Iₜ, Iₚ)
	n = size(W,2)
	edges = W * Iₚ # TF edges cannot be silent
	edges[.!silent_nodes(W, Iₜ),:] .= 0
	edges .!= 0
end
silent_edges(W, nₜ::Integer, nₚ::Integer) = silent_edges(W, Model.Iₜ(size(W,1),nₜ,nₚ), Model.Iₚ(size(W,1),nₜ,nₚ))
"Nodes where ongoing edges will be silent."
function silent_nodes(W, Iₜ)
	n = size(W,1)
	W′ = (W .!= 0)'
	(I - W′)^-1 * Iₜ * W′ * ones(n) .== 0
end
silent_nodes(W, nₜ::Integer, nₚ::Integer) = silent_nodes(W, Model.Iₜ(size(W,1),nₜ,nₚ))


"""
Proteins regulated only by phosphatases, which means the phosphate edges serve no purpose since the protein will always be unphosphorylated.
return: 1D bit vector, length size(Wₚ,1). true for each protein regulated by only phosphatases (and at least 1), false otherwise.
"""
function silent_phosphatases(Wₚ)
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
function correct_self_loops!(Wₜ, Wₚ)
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
function correct_silent_edges!(Wₜ, Wₚ)
	_, nₜ, nₚ = nₓnₜnₚ(Wₜ, Wₚ)
	W = Model._W(Wₜ, Wₚ)
	edges = silent_edges(W, nₜ, nₚ)
	check(edges, "silent edges")
	W[edges] .= 0
	Wₜ[:], Wₚ[:] = Model.WₜWₚ(W, nₜ, nₚ)
	return any(edges)
end
"""
Set (PK/PP) phosphate edges to zero if they will be silent in simulation and show a warning if any are found.
return: bool indicating if any were found.
"""
function correct_silent_phosphates!(Wₚ)
	phosphates = silent_phosphatases(Wₚ)
	check(phosphates, "silent phosphate regulations")
	Wₚ[phosphates,:] .= 0
	return any(phosphates)
end
"""
Correct self loops, silent edges, and silent phosphorylation regulations.
return: bool indicating if anything was corrected.
"""
function correct!(Wₜ, Wₚ)
	corrected = correct_self_loops!(Wₜ, Wₚ)
	while correct_silent_edges!(Wₜ, Wₚ) | correct_silent_phosphates!(Wₚ)
		corrected = true
	end
	corrected
end

end;
