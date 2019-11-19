#!/usr/bin/env julia
if !isdefined(Main, :Model) include("Model.jl") end
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end

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
Randomly split edges into sets P, T or X where X is nodes without outgoing edges, 
P has to have edges in B to nodes that is also regulated by T.
- nₚ: assign this many nodes to the set P.
"""
random_PTX(B::Matrix, nₚ::Integer) = random_PTX(B .!= 0, nₚ)
function random_PTX(B::BitMatrix, nₚ::Integer)
	n = size(B,1)
	# a node with zero outgoing edges belong to X
	X = vec(sum(B,dims=1) .== 0)
	# if a node is the only regulator of another, it will have to be in T,
	# since there would be no way for it to perform its regulation otherwise
	T = only_regulator(B)
	# find edges that are a "subset" as a PK edge: the target does not control more than the source, 
	# that is, bᵢⱼ is allowed if bₖⱼ >= bₖᵢ for any k.
	subset = subset_edges(B)
	# good candidates are all the nodes that has at least one "non-destructive" outgoing edge
	P_candidates = .!(X .| T) .& vec(sum(subset, dims=1) .!= 0)
	P = falses(n)
	# find nodes to add to P one at a time
	for i ∈ 1:nₚ
		# add one node to P at a time from the good candidates if there are any
		if any(P_candidates)
			P[rand((1:n)[P_candidates])] = true
		else  # otherwise add something random that changes B a bit
			remain = .!(P .|T .|X)
			if !any(remain)
				@warn("It was not possible to add as many nodes to P as wanted")
				break
			end
			P[rand((1:n)[remain])] = true
		end
		# T grows as well so we always make sure that a P does not end up being the only regulator of some node
		T[.!P] .|= only_regulator(B[:,.!P])
		P_candidates .&= .!(P .| T)  # remove anything in P or T from candidates
	end
	P, .!(P .|X), X
end

"""
Rewire P edges to random Ts that regulate their final target given as a B edge
"""
rewire(B::Matrix, P,T,X) = rewire(B .!= 0, P,T,X)
function rewire(B::BitMatrix, P,T,X)
	n = size(B,1)
	W = falses(n,n)
	W[:,T] .= B[:,T]  # T edges are the same
	P_num = (1:n)[P]  # logical to numerical index
	T_num = (1:n)[T]  # logical to numerical index
	for j ∈ P_num  # each P
		for i ∈ (1:n)[B[:,j]]  # each P b edge, numerical index
			# randomly select a tf to re-route through
			tf = rand(T_num[B[i,T]])
			W[tf,j] = true
			# check if we have fully described the intended effect (NO transcription cycles)
			if all(W * W[:,j] >= B[:,j]) break end
		end
	end
	W
end

"""
Find Ps that regulate "subsets" of what other Ps regulate and connect a random pair so that the one regulating less is the target.
Repeat until there is no more cascade edges to add.
"""
function add_cascades!(W, P)
	nₚₖ = sum(P)
	# which P can be a "subset" of another P?
	subset = subset_edges(W[:,P])
	while any(subset)
		# choose a random "subsetting" edge to add (linear index)
		idx = rand((1:nₚₖ^2)[vec(subset)])
		view(W,P,P)[idx] = true  # add P->P edge
		cartesian = CartesianIndices(subset)[idx]  # we need to know from and to which protein
		view(W,:,P)[:,cartesian[2]] .-= view(W,:,P)[:,cartesian[1]]  # remove edges that are now described indirectly through the new P->P edge
		# update
		subset = subset_edges(W[:,P])
	end
	W
end

"version where cascades are only added for half the places where it is possible"
function add_cascades50!(W, P)
	nₚₖ = sum(P)
	# which P can be a "subset" of another P?
	subset = subset_edges(W[:,P])
	n_subs = sum(subset)
	while sum(subset) > .5n_subs
		# choose a random "subsetting" edge to add (linear index)
		idx = rand((1:nₚₖ^2)[vec(subset)])
		view(W,P,P)[idx] = true  # add P->P edge
		cartesian = CartesianIndices(subset)[idx]  # we need to know from and to which protein
		view(W,:,P)[:,cartesian[2]] .-= view(W,:,P)[:,cartesian[1]]  # remove edges that are now described indirectly through the new P->P edge
		# update
		subset = subset_edges(W[:,P])
	end
	W
end

add_cascades_funs = Dict("100"=>add_cascades!, "050"=>add_cascades50!)

"""
Set a number of nodes ∈ P as a node ∈ PP, with the requirement that it cannot be the only phos regulator of any node.
Random extra kinase edges are added to balance out regulation.
- nₚₚ: number of phosphatases to create. 
"""
function phosphatases(W, P, nₚₚ::Integer)
	n = size(W,1)

	PP = tological(shuffle((1:n)[P])[1:nₚₚ], n)
	PK = P .& .!PP

	# find the nodes that has a PP regulator but no PK regulator
	only_has_PP = any(W[:,PP]; dims=2) .& .!any(W[:,PK]; dims=2)
	# add a single random PK edge onto each node that is currently only regulated by PP.
	# make sure it is not a self-loop that is added.
	for i ∈ (1:n)[vec(only_has_PP)]
		available_PKs = copy(PK)
		available_PKs[i] = false
		W[i,rand((1:n)[available_PKs])] = 1
	end
	
	# set the sign for PP
	W = 1W
	W[:,PP] .*= -1
	return W
end


"""
Set about half the TF edges as repressing.
- W: edges.
- T: logical index of which nodes are Ts.
"""
function repressors(W, T)
	W = 1W  # convert to int
	W[:,T] .*= rand([-1, 1], size(W))[:,T]
	W
end

"Sort an adjacency matrix so proteins are in order PK, PP, T, X."
function sort_PTX(W, P, T, X)
	n = size(W,1)
	order = [(1:n)[P]; (1:n)[T]; (1:n)[X]]
	reorder(W, order)
end

"Random Wₜ, Wₚ from B."
function random_W(B, nₚₖ::Integer, nₚₚ::Integer)
	P, T, X = random_PTX(B, nₚₖ + nₚₚ)
	W = rewire(B, P,T,X)
	add_cascades!(W, P)
	W = phosphatases(W, P, nₚₚ)
	W = repressors(W, T)
	W = sort_PTX(W, P, T, X)
	Model.WₜWₚ(W, sum(T), sum(P))
end

function random_W(B, nₚₖ::Integer, nₚₚ::Integer, fun)
	P, T, X = random_PTX(B, nₚₖ + nₚₚ)
	W = rewire(B, P,T,X)
	add_cascades_funs[fun](W, P)
	W = phosphatases(W, P, nₚₚ)
	W = repressors(W, T)
	W = sort_PTX(W, P, T, X)
	Model.WₜWₚ(W, sum(T), sum(P))
end

"Get a vector with -1 and 1 indicating which nodes ∈ P are ∈ PP and ∈ PK. 0 means it is in neither."
function PKPP(Wₚ::Matrix)
	PK = all(Wₚ .>= 0; dims=1)
	PP = all(Wₚ .<= 0; dims=1)
	(PK .& .!PP) - (PP .& .!PK)
end


end;

