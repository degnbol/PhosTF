#!/usr/bin/env julia
module FluxUtils
using LinearAlgebra
using Flux
import Formatting

export zerodiag
export random_weight

"""
Set Gaussian random values where σ is set by the number of values.
"""
random_weight(n::Integer, m::Integer) = randn(n, m) ./ √(n*m)

function Formatting.printfmt(prec::Integer, args::Real...; sep="\t")
    Formatting.printfmt(join(["{:.$(prec)f}" for _ in args], sep), args...)
end
Formatting.printfmtln(fe::AbstractString, args::Real...) = (Formatting.printfmt(fe, args...); println())
Formatting.printfmtln(prec::Integer, args::Real...; sep="\t") = (Formatting.printfmt(prec, args...; sep=sep); println())

end;
