#!/usr/bin/env julia
# First remove record of tested versions to do a clean install.
rm("Project.toml")
rm("Manifest.toml")

using Pkg

Pkg.add(url="https://github.com/mortenpi/ProjectX.jl.git")
Pkg.add("Revise")

Pkg.activate(".")

# precompile major packages to make sure they work before running jobs.
Pkg.add("DifferentialEquations")
using DifferentialEquations
Pkg.add("Flux")
using Flux

Pkg.add([
"ArgParse",
"BSON",
"CSV",
"Chain",
"Colors",
"DataFrames",
"Dates",
"DelimitedFiles",
"DiffEqCallbacks",
"Distributions",
"Fire",
"Formatting",
"Glob",
"HDF5",
"InlineStrings",
"IterTools",
"JLD",
"JSON3",
"LanguageServer",
"LinearAlgebra",
"Logging",
"NamedTupleTools",
"Random",
"SparseArrays",
"Statistics",
"StatsBase",
"Test",
"Plots",
"PlotlyJS",
])

Pkg.add(url="https://github.com/diegozea/ROC.jl")

# Using ProjectX we can have project local startup code auto loaded.
run(pipeline(`cat`; stdin="Project.toml.projectx", stdout="Project.toml", append=true))

