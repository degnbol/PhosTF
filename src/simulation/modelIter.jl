#!/usr/bin/env julia
include("ModelIteration.jl")
@use "utilities/ReadWrite"
using Fire

"""
Iterate model until convergence.
USE: ./modelIter.jl iter /path/to/WT.mat /path/to/WP.mat > outfile
"""
@main function iter(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = Main.ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

