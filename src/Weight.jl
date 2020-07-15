#!/usr/bin/env julia
isdefined(Main, :Model) || include("Model.jl")
isdefined(Main, :ArrayUtils) || include("utilities/ArrayUtils.jl")

"Construction of valid and random weight matrices, detection of silent edges and weight matrix correction."
module Weight
using LinearAlgebra
using SparseArrays
using ..Model, ..ArrayUtils
using Random


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
		E[V,1:nₚ] .|= W[V,1:nₚ] .!= 0   # mark KP edges as silent if they are onto a silent node
		V = vec(all((W .== 0) .| E, dims=1))  # recalculate silent nodes where the edges marked silent are removed from W
	end
	E
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
Set silent KP edges to zero and show a warning if any are found.
return: bool indicating if any were found.
"""
function correct_silent_edges!(W::AbstractMatrix, nₚ::Integer)
	E = silent_edges(W, nₚ)
	check(E, "silent edges")
	W[E] .= 0
	return any(E)
end
function correct_silent_edges!(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	_, nₜ, nₚ = Model.nₒnₜnₚ(Wₜ, Wₚ)
	W = Model._W(Wₜ, Wₚ)
	out = correct_silent_edges!(W, nₚ)
	Wₜ[:], Wₚ[:] = Model.WₜWₚ(W, nₜ, nₚ)
	return out
end

"""
Correct self loops, silent edges, and silent phosphorylation regulations.
return: bool indicating if anything was corrected.
"""
function correct!(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix)
	corrected = correct_self_loops!(Wₜ, Wₚ)
	while correct_silent_edges!(Wₜ, Wₚ) corrected = true end
	corrected
end
function correct!(W::AbstractMatrix, nₚ::Integer)
	corrected = correct_self_loops!(W)
	while correct_silent_edges!(W, nₚ) corrected = true end
	corrected
end

threshold!(W::AbstractMatrix, thres::AbstractFloat=0.001) = W[abs.(W) .< thres] .= 0



"""
Get an adjacency matrix where an edge indicates that target is a subset of source, 
in the sense of which outgoing edges each node has.
"""
function subset_edges(B::BitMatrix)
	out = B'B .>= sum(B, dims=1)'
	out[diagind(out)] .= 0  # no self-loops
	out[vec(sum(B,dims=1) .== 0),:] .= 0  # no subset defined for empty sets
	out
end

"Get nodes that only has a single regulator"
has_single_regulator(B) = vec(sum(B, dims=2) .== 1)
"Get the nodes which are the only regulator of some other node."
only_regulator(B::BitMatrix) = vec(sum(B[has_single_regulator(B),:], dims=1) .!= 0)

"""
Randomly split edges into sets KP, TF or O where O is nodes without outgoing edges, 
KP has to have edges in B to nodes that is also regulated by TF.
- nₚ: assign this many nodes to the set KP.
"""
random_KPTFO(B::Matrix, nₚ::Integer) = random_KPTFO(B .!= 0, nₚ)
function random_KPTFO(B::BitMatrix, nₚ::Integer)
	n = size(B,1)
	# a node with zero outgoing edges belong to O
    O = vec(sum(B,dims=1) .== 0)
	# if a node is the only regulator of another, it will have to be in TF,
	# since there would be no way for it to perform its regulation otherwise
	TF = only_regulator(B)
	# find edges that are a "subset" as a KP edge: the target does not control more than the source, 
	# that is, bᵢⱼ is allowed if bₖⱼ >= bₖᵢ for any k.
	subset = subset_edges(B)
	# good candidates are all the nodes that has at least one "non-destructive" outgoing edge
	KP_candidates = .!(O .| TF) .& vec(sum(subset, dims=1) .!= 0)
    
	KP = falses(n)
	# find nodes to add to KP one at a time
	for i ∈ 1:nₚ
		# add one node to KP at a time from the good candidates if there are any
		if any(KP_candidates)
			KP[rand((1:n)[KP_candidates])] = true
		else  # otherwise add something random that changes B a bit
			remain = .!(KP .|TF .|O)
			if !any(remain)
				@warn("It was not possible to add as many nodes to KP as wanted")
				break
			end
			KP[rand((1:n)[remain])] = true
		end
		# TF grows as well so we always make sure that a KP does not end up being the only regulator of some node
		TF[.!KP] .|= only_regulator(B[:,.!KP])
		KP_candidates .&= .!(KP .| TF)  # remove anything in KP or TF from candidates
	end
	KP, .!(KP .|O), O
end

"""
Rewire KP edges to random TFs that regulate their final target given as a B edge
"""
rewire(B::Matrix, KP, TF, O) = rewire(B .!= 0, KP, TF, O)
function rewire(B::BitMatrix, KP, TF, O)
	n = size(B,1)
	W = falses(n,n)
	W[:,TF] .= B[:,TF]  # TF edges are the same
	KP_num = (1:n)[KP]  # logical to numerical index
	TF_num = (1:n)[TF]  # logical to numerical index
	for j ∈ KP_num  # each KP
		for i ∈ (1:n)[B[:,j]]  # each KP b edge, numerical index
			# randomly select a TF to re-route through
			tf = rand(TF_num[B[i,TF]])
			W[tf,j] = true
			# check if we have fully described the intended effect (NO transcription cycles)
			all(W * W[:,j] >= B[:,j]) && break
		end
	end
	W
end

"""
Find KPs that regulate "subsets" of what other KPs regulate and connect a random pair so that the one regulating less is the target.
Repeat until there is no more cascade edges to add.
"""
function add_cascades!(W, KP)
	nₚ = sum(KP)
	# which KP can be a "subset" of another KP?
	subset = subset_edges(W[:,KP])
	while any(subset)
		# choose a random "subsetting" edge to add (linear index)
		idx = rand((1:nₚ^2)[vec(subset)])
		view(W,KP,KP)[idx] = true  # add KP->KP edge
		cartesian = CartesianIndices(subset)[idx]  # we need to know from and to which protein
		view(W,:,KP)[:,cartesian[2]] .-= view(W,:,KP)[:,cartesian[1]]  # remove edges that are now described indirectly through the new KP->KP edge
		# update
		subset = subset_edges(W[:,KP])
	end
	W
end

"""
Set about half the edges as repressing/deactivating.
- W: edges.
"""
function set_signs(W)
	W = 1W  # convert to int
	W .*= rand([-1, 1], size(W))
	W
end

"Sort an adjacency matrix so proteins are in order KP, TF, O."
function sort_KPTFO(W, KP, TF, O)
	n = size(W,1)
	order = [(1:n)[KP]; (1:n)[TF]; (1:n)[O]]
	reorder(W, order)
end

"Random Wₜ, Wₚ from B."
function random_W(B::AbstractMatrix, nₚ::Integer)
	KP, TF, O = random_KPTFO(B, nₚ)
	W = rewire(B, KP, TF, O)
	add_cascades!(W, KP)
	W = set_signs(W)
	W = sort_KPTFO(W, KP, TF, O)
	Model.WₜWₚ(W, sum(TF), sum(KP))
end


end;

