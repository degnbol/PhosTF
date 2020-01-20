#!/usr/bin/env julia
include("src/utilities/ReadWrite.jl")
include("src/utilities/CLI.jl")
include("src/utilities/General.jl")
include("src/ModelIteration.jl")
isdefined(Main, :Model) || include("src/Model.jl")
isdefined(Main, :ArrayUtils) || include("src/utilities/ArrayUtils.jl")


using Fire
using LinearAlgebra
using ..ReadWrite, ..ArrayUtils, ..General
using ..ModelIteration, ..Model
using ..CLI


"Iterate model until convergence."
@main function iter(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

"""
Get the total effects T from single KO experiments.
"""
@main function T(io=nothing, o=nothing)
	i, o = inout(io, o)
	X = loaddlm(i)
	savedlm(o, Model.X2T(X))
end

