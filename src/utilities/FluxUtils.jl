#!/usr/bin/env julia
module FluxUtils
using LinearAlgebra
using Flux
import Formatting

export zerodiag
export random_weight
export L1, L2

"""
Set Gaussian random values where σ is set by the number of values.
"""
random_weight(n::Integer, m::Integer) = randn(n, m) ./ √(n*m)
random_weight(ns::Tuple) = randn(ns...) ./ √prod(ns)

function Formatting.printfmt(prec::Integer, args::Real...; sep="\t")
    Formatting.printfmt(join(["{:.$(prec)f}" for _ in args], sep), args...)
end
Formatting.printfmtln(fe::AbstractString, args::Real...) = (Formatting.printfmt(fe, args...); println())
Formatting.printfmtln(prec::Integer, args::Real...; sep="\t") = (Formatting.printfmt(prec, args...; sep=sep); println())

# most efficient L1 and L2 regularization https://fluxml.ai/Flux.jl/stable/models/regularisation/
L1(x) = sum(abs, x)
L2(x) = sum(abs2, x)

end;
