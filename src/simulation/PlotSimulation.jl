#!/usr/bin/env julia

"PhosTF project plotting functions to visualize gene expression simulation etc."
module PlotSimulation
using Plots
using DifferentialEquations: ODESolution



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


"""
Plot simulation for multiple nodes.
- time: vector of time points for measurements.
- array: node values along axis 1, time along axis 2.
- labels for each column in array.
"""
function plot_simulation(time::Vector, values::Matrix, labels::Vector, styles::Vector, widths::Vector)
	p = plot()
	for i ∈ 1:length(labels)
		series = values[i,:]
		if all(series .== 0) continue end
		plot!(time, series, color="black", label=labels[i], style=styles[i], width=widths[i])
	end
	ylabel!("value/max")
	p
end


"""
Plot subplots each with their own set of curves all with a shared time axis. 
i-th subplot uses i-th element of values, labels, styles, widths and names.
"""
function plot_simulation(time::Vector, values::Vector{<:Matrix}, labels::Vector{<:Vector}, styles::Vector{<:Vector}, widths::Vector{<:Vector}, names::Vector{String})
	n_subplots = length(names)
	subplots = [plot_simulation(time, values[i], labels[i], styles[i], widths[i]) for i ∈ 1:n_subplots]
	plot(subplots..., layout=(n_subplots,1))
end


end;
