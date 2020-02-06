#!/usr/bin/env julia
include("src/utilities/ReadWrite.jl")
include("src/utilities/General.jl")
isdefined(Main, :Model) || include("src/Model.jl")

using Fire
using LinearAlgebra
using .ReadWrite, .General, .Model


function _load(Wₜ, Wₚ)
	Wₜ = loaddlm(abspath_(Wₜ), Float64)
	Wₚ = loaddlm(abspath_(Wₚ), Float64)
	Model._W(Wₜ, Wₚ), Model.constants(Wₜ, Wₚ)
end

@main function B(Wₜ, Wₚ, out="B.mat")
	savedlm(out, Model._B(_load(Wₜ, Wₚ)...))
end

@main function B_star(Wₜ, Wₚ, out="B_star.mat")
	savedlm(out, Model.B_star(_load(Wₜ, Wₚ)...))
end

