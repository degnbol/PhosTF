#!/usr/bin/env julia
module PlotTimeSeries
using DifferentialEquations: ODESolution
using DataFrames
using Plots
plotlyjs()


"""
Plot timeseries for multiple nodes.
Time along x-axis.
kwargs just there to (silently) throw away keyword args that aren't recognized.
"""
function plot_timeseries(df::DataFrame)
    if "subplot" ∉ names(df)
        subplots = [_plot_timeseries(df)]
    else
        subplots = [_plot_timeseries(df[df[!, "subplot"] .== v, :]) for v ∈ unique(df[!, "subplot"])]
    end
    xlabel!(subplots[end], "time [min]")
    plot(subplots..., layout=(length(subplots), 1))
end

function _plot_timeseries(df::DataFrame)
    p = plot()
    for i ∈ unique(df[!, "series"])
        series = df[df[!, "series"] .== i, :]
        # for a series all color, label, style, width should be the same so we use "first"
        kwargs = (Symbol(k)=>first(series[!, k]) for k ∈ ["color", "label", "style", "width"] if k ∈ names(series))
        plot!(series[!, "time"], series[!, "value"]; kwargs...)
	end
	ylabel!("value/max")
    ylims!(0, 1)
	p
end

end;
