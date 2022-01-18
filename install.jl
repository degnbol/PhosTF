#!/usr/bin/env julia
using Pkg
if ARGS[1] == "main"
    println("installs julia dependencies to main environment.")
    Pkg.add("TOML")
    using TOML
    dep_names = TOML.parsefile("Project.toml")["deps"] |> keys .|> String
    Pkg.add(dep_names)
elseif ARGS[1] == "env"
    println("installs julia dependencies to an isolated environment.")
    Pkg.activate(".")
    Pkg.instantiate()
else
    error("Supply single argument \"main\" or \"env\".")
end 
