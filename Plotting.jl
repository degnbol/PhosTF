#!/usr/bin/env julia

"PKTFX project plotting functions to visualize gene expression simulation etc."
module Plotting
using Plots
using ColorSchemes
using DifferentialEquations: ODESolution


color_scheme = ColorSchemes.phase


function save(o, p)
	if o == stdout display(p)
	else
		ext = splitext(o)[2]
		if ext == ".png" png(o)
		elseif ext == ".pdf" savefig(p, o)
		else error("Unrecognized plot file format.") end
	end
end


plot_r(solution::ODESolution) = plot(solution.t, solution[:,1,:]', leg=false)
plot_p(solution::ODESolution) = plot(solution.t, solution[:,2,:]', leg=false)
plot_ϕ(solution::ODESolution) = plot(solution.t, solution[:,3,:]', leg=false)
plot_ψ(solution::ODESolution) = plot(solution.t, solution[:,2,:]'-solution[:,3,:]', leg=false)
plot_ϕ_frac(solution::ODESolution) = plot(solution.t, solution[:,3,:]'./solution[:,2,:]', leg=false)


plot_simulation(time, array::Matrix, color) = plot(time, array', color=color, leg=false)
plot_simulation!(time, array::Matrix, color) = plot!(time, array', color=color, leg=false)
function plot_simulation(time, array::Vector{Matrix{Float64}})
	n = length(array)
	p = plot_simulation(time, array[1], get(color_scheme, 1/(n+1)))
	for i in 2:n
		p = plot_simulation!(time, array[i], get(color_scheme, i/(n+1)))
	end
	p
end


end;