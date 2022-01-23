#!/usr/bin/env julia
"Plotting time series."
module PlotTimeSeries
using Plots
using DifferentialEquations: ODESolution

"""
Plot timeseries for multiple nodes.
- time: vector of time points for measurements.
- values: node values along axis 1, time along axis 2.
- labels for each column in array.
"""
function plot_timeseries(time::Vector, values::Matrix, labels::Vector, styles::Vector, widths::Vector; colors=nothing)
	p = plot()
    if colors === nothing colors = ["black" for l in labels] end
    for (i, col) ∈ enumerate(colors)
		series = values[i, :]
		any(series .!= 0) || continue
		plot!(time, series, color=col, label=labels[i], style=styles[i], width=widths[i])
	end
	ylabel!("value/max")
	p
end


"""
Plot subplots each with their own set of curves all with a shared time axis. 
i-th subplot uses i-th element of values, labels, styles, widths and names.
- values: Vector{Matrix{Float}} but using that in the declaration seems to be too specific.
"""
function plot_timeseries(time::Vector, values::Vector, labels::Vector{Vector{String}}, styles::Vector{Vector{Symbol}}, widths::Vector{Vector{Int}}, names::Vector{String}; colors=nothing)
#= function plot_timeseries(time::Vector, values::Vector, labels::Vector, styles::Vector, widths::Vector, names::Vector; colors=nothing) =#
	n_subplots = length(names)
	subplots = [plot_timeseries(time, values[i], labels[i], styles[i], widths[i]; colors=colors) for i ∈ 1:n_subplots]
	plot(subplots..., layout=(n_subplots, 1))
end

end;
