#!/usr/bin/env julia
include("../utilities/ReadWrite.jl")
include("../utilities/CLI.jl")
include("ModelIteration.jl")
isdefined(Main, :Model) || include("../inference/Model.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")


using Fire
using LinearAlgebra
using ..ReadWrite, ..ArrayUtils
using ..ModelIteration, ..Model
using ..CLI


"""
Get the total effects T from single KO experiments.
Used in comparison with LLC_X2T.R in etc/LLC to validate correct calculation of T.
USE: ./X2T.jl T infile outfile
"""
@main function T(io=nothing, o=nothing)
	i, o = inout(io, o)
	X = loaddlm(i)
	savedlm(o, Model.X2T(X))
end

