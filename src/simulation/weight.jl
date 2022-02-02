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
- nₚ: how many of the nodes should be attempted to be assigned as P rather than T?
- pretty: set to true to get outputs written with '.', '+', '-' instead of 0, 1, -1.
"""
@main function random(B::String, nₚ::Integer; WT::String="WT.adj", WP::String="WP.adj", pretty::Bool=false)
	Wₜ, Wₚ = WeightConstruction.random_WₜWₚ(loaddlm(B), nₚ)
    if pretty
        Wₜ = ReadWrite.pretty_matrix(Wₜ)
        Wₚ = ReadWrite.pretty_matrix(Wₚ)
    end
	# the nodes gets sorted TF, KP, O so we can add some row and column names
    nₜ, nₚ, nₒ = Model.nₜnₚnₒ(Wₜ, Wₚ)
    names = vcat(["T$i" for i in 1:nₜ], ["P$i" for i in 1:nₚ], ["O$i" for i in 1:nₒ])
    savedlm(WT, Wₜ; colnames=names[1:nₜ], rownames=names)
    savedlm(WP, Wₚ; colnames=names[nₜ+1:nₜ+nₚ], rownames=names[1:nₜ+nₚ])
end

"""
This version of the function is run when no arguments are supplied.
"""
@main function correct(Wₜfname::String="WT.adj", Wₚfname::String="WP.adj", ot="WT_cor.adj", op="WP_cor.adj", save::Bool=false)
    Wₜ = loadmat(Wₜfname; header=true)
    Wₚ = loadmat(Wₚfname; header=true)
    
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
- header: does the matrix have a header (and potentially row names)
- thres: optional threshold
"""
@main function thres(io=nothing, o=nothing; header::Bool=false, thres=0.001)
	i, o = inout(io, o)
	mat = loaddlm(i, Float64; header=header)
	WeightConstruction.threshold!(mat, thres)
	savedlm(o, mat)
end

