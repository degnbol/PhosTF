#!/usr/bin/env julia
isdefined(Main, :ReadWrite) || include("../utilities/ReadWrite.jl")
isdefined(Main, :CLI) || include("../utilities/CLI.jl")
isdefined(Main, :General) || include("../utilities/General.jl")
isdefined(Main, :WeightConstruction) || include("WeightConstruction.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")

using Fire
using .ReadWrite, .ArrayUtils, .General, .WeightConstruction, .CLI

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"

"""
Create random W from a adjacency matrix containing B.
- B: filename for B matrix, usually goldstandard matrix from DREAM4.
- nₚ: how many of the nodes should be attempted to be assigned as KP rather than TF?
"""
@main function random(B::String, nₚ::Integer; WT_fname::String=default_Wₜ, WP_fname::String=default_Wₚ)
	Wₜ, Wₚ = WeightConstruction.random_W(loaddlm(B), nₚ)
	savedlm(WT_fname, Wₜ)
	savedlm(WP_fname, Wₚ)
end

"""
- np: mandatory arg to set nₚ
- save: should we save files if no corrections are made?
"""
@main function correct(io=nothing, o=nothing; np=nothing, save::Bool=false)
	if np === nothing @error("supply --np"); return end
	i, o = inout(io, o)
	W = loaddlm(i)
	if WeightConstruction.correct!(W, np)
		@info("Corrections made.")
		savedlm(o, W)
	else
		@info("NO corrections made.")
		save && savedlm(o, W)
	end
end
"""
This version of the function is run when no arguments are supplied.
"""
@main function correct(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; ot="WT_cor.mat", op="WP_cor.mat", save::Bool=false)
	Wₜ, Wₚ = loadmat(Wₜ_fname), loadmat(Wₚ_fname)

	if WeightConstruction.correct!(Wₜ, Wₚ)
		@info("Corrections made.")
		savedlm(ot, Wₜ)
		savedlm(op, Wₚ)
	else
		@info("NO corrections made.")
		if save
			savedlm(ot, Wₜ)
			savedlm(op, Wₚ)
		end
	end
end

"""
Remove edges less than a given threshold.
- thres: optional threshold
"""
@main function thres(io=nothing, o=nothing; thres=0.001)
	i, o = inout(io, o)
	mat = loaddlm(i, Float64)
	WeightConstruction.threshold!(mat, thres)
	savedlm(o, mat)
end


"""
Swap order of KP and TF in matrix (swap between KP-TF-V and TF-KP-V).
- n1: number of nodes in first group, which will be moved to become the second.
- n2: number of nodes in the second group, which will be move to become the first.
"""
@main function swapPT(io=nothing, o=nothing; n1=nothing, n2=nothing)
	i, o = inout(io, o)
	if n1 === nothing || n2 === nothing @error("Provide --n1 and --n2.") end
	mat = loaddlm(i, Float64)
	mat = reorder(mat, [n1+1:n1+n2;1:n1;n1+n2+1:maximum(size(mat))])
	savedlm(o, mat)
end

