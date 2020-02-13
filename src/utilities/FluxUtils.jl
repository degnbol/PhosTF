#!/usr/bin/env julia
module FluxUtils
using LinearAlgebra
using Flux
using Flux: Tracker
import Formatting

export zerodiag
export random_weight

"""
Set Gaussian random values where σ is set by the number of values.
"""
random_weight(n::Integer, m::Integer) = randn(n, m) ./ √(n*m)

Formatting.printfmt(fe::AbstractString, args::Tracker.TrackedReal...) = Formatting.printfmt(fe, [Tracker.data(v) for v in args]...)
function Formatting.printfmt(prec::Integer, args::Tracker.TrackedReal...; sep="\t")
    Formatting.printfmt(join(["{:.$(prec)f}" for _ in args], sep), args...)
end
Formatting.printfmtln(fe::AbstractString, args::Tracker.TrackedReal...) = (Formatting.printfmt(fe, args...); println())
Formatting.printfmtln(prec::Integer, args::Tracker.TrackedReal...; sep="\t") = (Formatting.printfmt(prec, args...; sep=sep); println())

end;