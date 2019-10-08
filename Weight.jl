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


function correct(Wₜ, Wₚ)
	nₓ, nₜ, nₚ = nₓnₜnₚ(Wₜ, Wₚ)
	self_loops_Wₜ, self_loops_Wₚ = self_loops(Wₜ), self_loops(Wₚ)
	check(self_loops_Wₜ, )
	check(self_loops_Wₚ)
end

end;
