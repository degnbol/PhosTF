#!/usr/bin/env julia
# First remove record of tested versions to do a clean install.
rm("Project.toml")
rm("Manifest.toml")

using Pkg
Pkg.activate(".")

Pkg.add("DifferentialEquations")
Pkg.add("Flux")

# precompile major packages to make sure they work before running jobs.
using Flux, DifferentialEquations

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
"Revise",
"SparseArrays",
"Statistics",
"StatsBase",
"Test",
])

Pkg.add(url="https://github.com/diegozea/ROC.jl")
Pkg.add(url="https://github.com/mortenpi/ProjectX.jl.git")

# Using ProjectX we can have project local startup code auto loaded.
run(`cat Project.toml.projectx >> ./Project.toml`)

