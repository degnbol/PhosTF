#!/usr/bin/env julia
# installs julia dependencies.
using Pkg
Pkg.activate("src")
Pkg.instantiate()
