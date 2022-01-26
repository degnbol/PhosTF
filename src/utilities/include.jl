#!/usr/bin/env julia
ROOT = readchomp(`git root`) * '/'
SRC = "$(ROOT)src/"
# change this flag in interactive mode when debugging to apply changes.
OVERWRITE = false

"""
Load a file into Main scope only providing path relative to root folder and no extension.
- path: e.g. "src/simulation/simulate"
"""
incl(path::String) = Base.include(Main, "$ROOT$path.jl")

"""
Load SOURCE file into Main scope (not global module scope) IF they aren't already defined.
paths: strings in the form "utilities/ArrayUtils" ... etc. i.e. relative to source folder.
"""
macro src(path::String)
    sym = path |> basename |> Symbol
    (!Main.OVERWRITE && isdefined(Main, sym)) || Base.include(Main, "$SRC$path.jl");
end
"Same as @incl but using is also called on the file (so it should contain a module)."
macro use(path::String)
    sym = path |> basename |> Symbol
    (!Main.OVERWRITE && isdefined(Main, sym)) || Base.include(Main, "$SRC$path.jl");
    :(using Main.$sym)
end

