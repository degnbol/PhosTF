#!/usr/bin/env julia
"Construction of valid and random weight matrices, detection of silent edges and weight matrix correction."
module WeightConstruction
using LinearAlgebra
using SparseArrays
Main.@use "inference/Model"
Main.@use "utilities/ArrayUtils"
using Random
using Chain: @chain
using DataFrames

"""
Get Wₜ, Wₚ from square W with nodes sorted in order TF, KP, O.
"""
WₜWₚ(W::AbstractMatrix, nₜ::Integer, nₚ::Integer) = W[:, 1:nₜ], W[1:nₜ+nₚ, nₜ+1:nₜ+nₚ]

self_loops(mat::AbstractMatrix, k::Integer=0) = (diag(mat, k) .!= 0)

"""
Mask of edges that leads to deadends, which are phosphorylation regulation leading to nodes that regulate nothing.
"""
function silent_edges(W::AbstractMatrix, nₚ::Integer)
	n = size(W, 1)
	O = vec(all(W .== 0, dims=1))  # nodes not regulating anything
	E = falses(n, n)
	E_last = trues(n,n)  # just any value that is different from E
	while any(E .!= E_last)
		E_last .= E
		E[O, 1:nₚ] .|= W[O, 1:nₚ] .!= 0   # mark KP edges as silent if they are onto a silent node
		O = vec(all((W .== 0) .| E, dims=1))  # recalculate silent nodes where the edges marked silent are removed from W
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
	nₚ = size(Wₚ, 2)
	self_loops_Wₜ, self_loops_Wₚ = self_loops(Wₜ, -nₚ), self_loops(Wₚ)
	check(self_loops_Wₜ, "self loops in Wₜ")
	check(self_loops_Wₚ, "self loops in Wₚ")
	Wₜ[diagind(Wₜ, -nₚ)[self_loops_Wₜ]] .= 0
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
    nₜ, nₚ = size(Wₜ, 2), size(Wₚ, 2)
	W = _W(Wₜ, Wₚ)
	out = correct_silent_edges!(W, nₚ)
	Wₜ[:], Wₚ[:] = WₜWₚ(W, nₜ, nₚ)
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
function correct!(Wₜ::DataFrame, Wₚ::DataFrame)
    @assert eltype(Wₜ[!,1]) <: AbstractString
    @assert eltype(Wₚ[!,1]) <: AbstractString
    R_names = [names(Wₜ)[2:end]; names(Wₚ)[2:end]]
    @assert all(R_names .== Wₚ[:,1] .== Wₜ[1:length(R_names),1])
    Wₜmat = Matrix(Wₜ[:, 2:end])
    Wₚmat = Matrix(Wₚ[:, 2:end])
    corrected = correct!(Wₜmat, Wₚmat)
    Wₜ[:, 2:end] .= Wₜmat
    Wₚ[:, 2:end] .= Wₚmat
    corrected
end

threshold!(W::AbstractMatrix, thres::AbstractFloat=0.001,) = W[abs.(W) .< thres] .= 0
threshold!(W::DataFrame, thres::AbstractFloat=0.001,) = begin
    for col ∈ eachcol(W)
        if eltype(col) <: Real
            col[col .< thres] .= 0
        end
    end
end


"""
Get an adjacency matrix where an entry indicates that the target (row) regulates a subset of what source (column) regulates.
Only direct regulation is considered, i.e. single edge.
Works for both BitMatrix and other Real Matrices, where all signs will have to match for the latter, for an element to be set to true.
return: square BitMatrix with both dimensions of length = size(adj,2), i.e. an entry from any source to any source
"""
function subset_adj(adj::AbstractMatrix)::BitMatrix
    nOut = sum(abs.(adj), dims=1)
    # for each pairwise combination of nodes, get the product of their outgoing edges with adj'adj.
    # If that is greater or equal to the total number of ougoing edges for one of the nodes it regulates a subset of the other.
    # It also needs to be true that nOut is larger than 0, i.e. nodes that regulates nothing is not regulating a subset of others.
	S = adj'adj .>= nOut' .> 0
	S[diagind(S)] .= 0  # no self-loops
	S
end

"Get nodes that only has a single regulator"
has_single_regulator(B) = vec(sum(B, dims=2) .== 1)
"Get the nodes which are the only regulator of some other node."
only_regulator(B::BitMatrix) = vec(sum(B[has_single_regulator(B), :], dims=1) .!= 0)

"""
Randomly split edges into sets KP, TF or O where O is nodes without outgoing edges, 
KP has to have edges in B to nodes that is also regulated by TF.
- nₚ: assign this many nodes to the set KP, or at least try to get this many valid nodes.
"""
random_TFKPO(B::Matrix, nₚ::Integer) = random_TFKPO(B .!= 0, nₚ)
function random_TFKPO(B::BitMatrix, nₚ::Integer)::Tuple{BitVector,BitVector,BitVector}
	n = size(B, 1)
	# a node with zero outgoing edges belong to O
    O = vec(sum(B, dims=1) .== 0)
	# if a node is the only regulator of another, it will have to be in TF,
	# since there would be no way for it to perform its regulation otherwise
	TF = only_regulator(B)
	# find edges that are a "subset" as a KP edge: the target does not control more than the source, 
	# that is, bᵢⱼ is allowed if bₖⱼ >= bₖᵢ for any k.
	B_sub = subset_adj(B)
	# good candidates are all the nodes that has at least one "non-destructive" outgoing edge
    KP_candidates = .!(TF .| O) .& vec(sum(B_sub, dims=1) .!= 0)
    
	KP = falses(n)
	# find nodes to add to KP one at a time
	for i ∈ 1:nₚ
		# add one node to KP at a time from the good candidates if there are any
		if any(KP_candidates)
			KP[rand((1:n)[KP_candidates])] = true
		else  # otherwise add something random that changes B a bit
			remain = .!(KP .| TF .| O)
			if !any(remain)
				@warn("It was not possible to add as many nodes to KP as wanted")
				break
			end
			KP[rand((1:n)[remain])] = true
		end
		# TF grows as well so we always make sure that a KP does not end up being the only regulator of some node
		TF[.!KP] .|= only_regulator(B[:, .!KP])
		KP_candidates .&= .!(TF .| KP)  # remove anything in TF or KP from candidates
	end
	.!(KP .|O), KP, O
end

"""
Rewire KP edges to random TFs that regulate their final target given as a B edge
"""
function redirect_KP_thru_TF(B::AbstractMatrix, TF::BitVector, KP::BitVector, O::BitVector)
	n = size(B, 1)
    W = copy(B) # it will have same type as B, e.g. BitMatrix or Matrix{Int}
    # TF->V edges are unchanged
	W[:, KP] .= 0
    for kp ∈ (1:n)[KP]  # each KP
		for target ∈ (1:n)[B[:, kp] .!= 0]  # each KP regulatory target. Numerical index
			# Find the best TF to reroute through. First find all TFs, where TF->target
            TFs_reg_target = (1:n)[(B[target, :] .!= 0) .& TF]
			# The weight that the needs to be used to result in the same product of sign through KP->TF->target as it was through KP->target.
            # In case of BitMatrix this will always just be true * true == true
            ws = B[target, kp] .* B[target, TFs_reg_target]
            # count how many values does not match between the KP and the TFs being considered.
            imperfections = sum(B[:, TFs_reg_target] .* ws' .!= B[:, kp], dims=1) |> vec
			# Then select the one that closest match the KP in its intended regulatory effect
            W[TFs_reg_target[argmin(imperfections)], kp] = ws[argmin(imperfections)]
			# check if we have fully described the intended effect (NO transcription cycles)
			all(W * W[:, kp] .!= 0 .>= B[:, kp]) && break
		end
	end
	W
end

"""
Find KPs that regulate subsets of the set of nodes that other KPs regulate. 
Then connect a random pair so that the one regulating less becomes the target of the other.
Repeat until there is no more cascade edges to add.
"""
function redirect_KP_thru_KP(W::AbstractMatrix, KP::BitVector)
    W = copy(W)
	nₚ = sum(KP)
	# which KP regulates a subset of the nodes that another KP regulates?
	Wₚsub = subset_adj(W[:, KP])
	while any(Wₚsub)
		# choose a random "subsetting" edge to add (converting to 1D index, aka. linear index)
		idx = rand((1:nₚ^2)[vec(Wₚsub)])
		view(W, KP, KP)[idx] = true  # add KP->KP edge
		cartesian = CartesianIndices(Wₚsub)[idx]  # we need to know from and to which protein
		view(W, :, KP)[:, cartesian[2]] .-= view(W, :, KP)[:, cartesian[1]]  # remove edges that are now described indirectly through the new KP->KP edge
		# update
		Wₚsub = subset_adj(W[:, KP])
	end
	W
end

"""
Set about half the edges as repressing/deactivating.
- W: TF and KP edges, unordered proteins are allowed.
"""
function set_signs(W::BitMatrix)
	W = 1W  # convert to int
	W .*= rand([-1, 1], size(W))
	W
end
set_signs(W::Matrix) = W

"Sort an adjacency matrix so proteins are in order TF, KP, O."
function sort_TFKPO(W::AbstractMatrix, TF::BitVector, KP::BitVector, O::BitVector)
	n = size(W, 1)
	order = [(1:n)[TF]; (1:n)[KP]; (1:n)[O]]
	reorder(W, order)
end

"Random Wₜ, Wₚ from B. Main function in this file."
function random_WₜWₚ(B::AbstractMatrix, nₚ::Integer)
	TF, KP, O = random_TFKPO(B, nₚ)
    @info "Found $(sum(O)) non-regulators assigned to O, $(sum(KP)) that could be assigned to KP, and the remaining $(sum(TF)) were assigned to TF."
	@chain B begin
	    redirect_KP_thru_TF(TF, KP, O)
	    redirect_KP_thru_KP(KP)
        set_signs
        sort_TFKPO(TF, KP, O)
        WₜWₚ(sum(TF), sum(KP))
    end
end


end;

