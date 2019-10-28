#!/usr/bin/env julia
include("Model.jl")

"Iterating the equations from the model intended for inference until convergence in order to generate naive data."
module ModelIteration
using Distributions: Normal, mean
using ..Model: Constants, nₓnₜnₚ, _W, _B

export converge

default_tolerance = 1e-7
default_max_iterations = 1000

function random_C(U, dist=Normal(-4))
	rand(dist, size(U)) .* (1 .- U)
end


"""
For iterating the model.
This is only for sanity checks and debugging since it imprecisely assumes equilibrium eq. B will work for each iteration, and naturally equilibrium assumption does not hold during each simulation step.
"""
step(X::Matrix, W, cs::Constants, C::Matrix) = _B(cs,W) * X .* cs.U .+ C
step(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = (_B(cs,W) * X + E) .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix) = X .= _B(cs,W) * X .* cs.U .+ C
step!(X::Matrix, W, cs::Constants, C::Matrix, E::Matrix) = X .= (_B(cs,W) * X + E) .* cs.U .+ C


function converge(W, constants::Constants, C=random_C(constants.U), X₀=C; tolerance=default_tolerance, max_iterations=default_max_iterations)
	Xₜ = Xₜ₋₁ = X₀
	difference = nothing
	for it in 1:max_iterations
		Xₜ = step(Xₜ₋₁, W, constants, C)
		difference = mean(abs.(Xₜ .- Xₜ₋₁))
		if difference < tolerance return Xₜ end
		Xₜ₋₁ = Xₜ
	end
	@error("Convergence was not reached ($difference !< $tolerance)")
	Xₜ
end
function converge(Wₜ, Wₚ, C=nothing, X₀=nothing; tolerance=default_tolerance, max_iterations=default_max_iterations)
	nₓ,nₜ,nₚ = nₓnₜnₚ(Wₜ, Wₚ)
	K = C === nothing ? nₜ+nₚ : size(C,2)
	constants = Constants(nₓ+nₜ+nₚ, nₜ, nₚ, K)
	if C === nothing C = random_C(constants.U) end
	if X₀ === nothing X₀ = C end
	converge(_W(Wₜ, Wₚ), constants, C, X₀, tolerance=tolerance, max_iterations=max_iterations)
end

end;