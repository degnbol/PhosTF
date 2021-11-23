#!/usr/bin/env julia
include("../utilities/ReadWrite.jl")
include("../utilities/CLI.jl")
include("../utilities/General.jl")
include("ModelIteration.jl")
isdefined(Main, :Model) || include("../inference/Model.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")


using Fire
using LinearAlgebra
using ..ReadWrite, ..ArrayUtils, ..General
using ..ModelIteration, ..Model
using ..CLI

"""
Iterate model until convergence.
USE: ./modelIter.jl iter /path/to/WT.mat /path/to/WP.mat > outfile
"""
@main function iter(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

