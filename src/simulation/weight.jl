#!/usr/bin/env julia
@use "simulation/WeightConstruction"
@src "inference/Model"
@use "utilities/ReadWrite"
@use "utilities/ArrayUtils"
@use "utilities/CLI"
using Fire

"""
Create random W from a adjacency matrix containing B.
- B: filename for B matrix, usually goldstandard matrix from DREAM4.
- nₚ: how many of the nodes should be attempted to be assigned as KP rather than TF?
- pretty: set to true to get outputs written with '.', '+', '-' instead of 0, 1, -1.
"""
@main function random(B::String, nₚ::Integer; WT::String="WT.mat", WP::String="WP.mat", pretty::Bool=false)
	Wₜ, Wₚ = WeightConstruction.random_WₜWₚ(loaddlm(B), nₚ)
    if pretty
        Wₜ = ReadWrite.pretty_matrix(Wₜ)
        Wₚ = ReadWrite.pretty_matrix(Wₚ)
    end
	# the nodes gets sorted TF, KP, O so we can add some row and column names
    nₜ, nₚ, nₒ = Model.nₜnₚnₒ(Wₜ, Wₚ)
    names = vcat(["TF$i" for i in 1:nₜ], ["KP$i" for i in 1:nₚ], ["O$i" for i in 1:nₒ])
    savedlm(WT, Wₜ; colnames=names[1:nₜ], rownames=names)
    savedlm(WP, Wₚ; colnames=names[nₜ+1:nₜ+nₚ], rownames=names[1:nₜ+nₚ])
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
@main function correct(Wₜfname::String="WT.mat", Wₚfname::String="WP.mat", op="WP_cor.mat", save::Bool=false)
	Wₜ, Wₚ = loadmat(Wₜfname), loadmat(Wₚfname)
    
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

