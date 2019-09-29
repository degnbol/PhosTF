#!/usr/bin/env julia
module ModelIteration
using Distributions: Normal, mean
include("Model.jl")

export converge

default_tolerance = 1e-7
default_max_iterations = 1000

function random_C(U, dist=Normal(-4))
	rand(dist, size(U)) .* (1 .- U)
end

function converge(W, constants::Model.Constants, C=random_C(constants.U), X₀=C; tolerance=default_tolerance, max_iterations=default_max_iterations)
	Xₜ = Xₜ₋₁ = X₀
	difference = nothing
	for it in 1:max_iterations
		Xₜ = Model.step(Xₜ₋₁, W, constants, C)
		difference = mean(abs.(Xₜ .- Xₜ₋₁))
		if difference < tolerance return Xₜ end
		Xₜ₋₁ = Xₜ
	end
	@error("Convergence was not reached ($difference !< $tolerance)")
	Xₜ
end
function converge(Wₜ, Wₚ, C=nothing, X₀=nothing; tolerance=default_tolerance, max_iterations=default_max_iterations)
	(n, nₜ), nₚ = size(Wₜ), size(Wₚ,2)
	K = C == nothing ? nₜ+nₚ : size(C,2)
	constants = Model.Constants(n, nₜ, nₚ, K)
	if C == nothing C = random_C(constants.U) end
	if X₀ == nothing X₀ = C end
	converge(Model._W(Wₜ, Wₚ), constants, C, X₀, tolerance=tolerance, max_iterations=max_iterations)
end

end;