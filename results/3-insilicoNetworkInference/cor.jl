#!/usr/bin/env julia
using DataFrames, CSV
using Statistics: cor, mean
using Glob
SRC = readchomp(`git root`) * "/src/"
include(SRC * "utilities/ReadWrite.jl")

readvecs(fname1, fname2) = begin
    df1 = CSV.read(fname1, DataFrame; delim=' ')
    df2 = CSV.read(fname2, DataFrame; delim=' ')
    @assert all(names(df1) == names(df2))
    @assert all(df1[!, 1] .== df2[!, 1])
    vec1 = ReadWrite.parse_matrix(Matrix(df1[!, 2:end])) |> vec
    vec2 = ReadWrite.parse_matrix(Matrix(df2[!, 2:end])) |> vec
    vec1, vec2
end
corr(fname1::AbstractString, fname2::AbstractString) = cor(readvecs(fname1, fname2)...)
avg(fnames1, dir) = [corr(f, dir * split(f, '/')[end]) for f in fnames1] |> mean


avg(Glob.glob("../*-insilicoNetworkSimulation/adjacencies/WT*.adj"), "inferredWeights/") |> println
avg(Glob.glob("../*-insilicoNetworkSimulation/adjacencies/WP*.adj"), "inferredWeights/") |> println


