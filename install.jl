#!/usr/bin/env julia
using Pkg
if length(ARGS) == 0
    println("installs julia dependencies to main environment.")
    Pkg.add("TOML")
    using TOML
    dep_names = TOML.parsefile("src/Project.toml")["deps"] |> keys .|> String
    Pkg.add(dep_names)
elseif ARGS[1] == "env"
    println("installs julia dependencies to an isolated environment.")
    Pkg.activate("src")
    Pkg.instantiate()
else
    error("Unknown option.")
end 
