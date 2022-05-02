#!/usr/bin/env julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# precompile major packages to make sure they work before running jobs.
using Flux, DifferentialEquations

