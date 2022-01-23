#!/usr/bin/env julia
@src "inference/Model"
@use "utilities/ReadWrite"
@src "utilities/CLI"
using Fire

"""
Get the total effects T from single KO experiments.
Used in comparison with LLC_X2T.R in etc/LLC to validate correct calculation of T.
USE: ./X2T.jl T infile outfile
"""
@main function T(io=nothing, o=nothing)
	i, o = CLI.inout(io, o)
	X = loaddlm(i)
	savedlm(o, Model.X2T(X))
end

