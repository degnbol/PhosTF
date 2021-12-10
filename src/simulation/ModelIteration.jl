#!/usr/bin/env julia
include("../inference/Model.jl")

"Iterating the equations from the model intended for inference until convergence in order to generate naive data."
module ModelIteration
using Distributions: Normal, mean
using ..Model

export converge

default_tolerance = 1e-7
default_max_iterations = 1000

function random_C(U, dist=Normal(-4))
    J = 1 .- U
	rand(dist, size(U)) .* J
end


"""
For iterating the model.
This is only for sanity checks and debugging since it imprecisely assumes equilibrium eq. B will work for each iteration, and naturally equilibrium assumption does not hold during each simulation step.
"""
step(mdl, X::Matrix, C::Matrix) = mdl(X) .* mdl.U .+ C
step(mdl, X::Matrix, C::Matrix, E::Matrix) = (Model._B(mdl) * X + E) .* mdl.U .+ C
step!(mdl, X::Matrix, C::Matrix) = X .= Model._B(mdl) * X .* mdl.U .+ C
step!(mdl, X::Matrix, C::Matrix, E::Matrix) = X .= (Model._B(mdl) * X + E) .* mdl.U .+ C


function converge(mdl, C=random_C(mdl.U), X₀=C; tolerance::Float64=default_tolerance, max_iterations::Integer=default_max_iterations)
	Xₜ = Xₜ₋₁ = X₀
	difference = nothing
	for it in 1:max_iterations
		Xₜ = step(mdl, Xₜ₋₁, C)
		difference = mean(abs.(Xₜ .- Xₜ₋₁))
		if difference < tolerance return Xₜ end
		Xₜ₋₁ = Xₜ
	end
	@error("Convergence was not reached ($difference !< $tolerance)")
	Xₜ
end
function converge(Wₜ::AbstractMatrix, Wₚ::AbstractMatrix, C=nothing, X₀=nothing; tolerance::Float64=default_tolerance, max_iterations::Integer=default_max_iterations)
	nₜ, nₚ, nₒ = Model.nₜnₚnₒ(Wₜ, Wₚ)
	K = C === nothing ? nₜ+nₚ : size(C, 2)
	mdl = Model.get_model(nₜ, nₚ, nₒ, K)
	if C === nothing C = random_C(mdl.U) end
	if X₀ === nothing X₀ = C end
	converge(mdl, C, X₀, tolerance=tolerance, max_iterations=max_iterations)
end

end;
