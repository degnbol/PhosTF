#!/usr/bin/env julia
include("GeneRegulation.jl")
"""
Defining and solving ODEs to the point of having the resulting simulated logFC values.
Used for simulation of gene expression levels.
"""
module ODEs
using DifferentialEquations: ODEProblem, solve, ODESolution
using DiffEqCallbacks: TerminateSteadyState
using Distributions: Uniform
using LinearAlgebra: I
using ..GeneRegulation

export simulate, steady_state
export logFC

""" ϕ₀ == 0 for non-regulating proteins in order to have vectors match in length """
get_u₀(net::Network) = [net.r₀ net.p₀ [net.ϕ₀; zeros(net.n - (net.nₚ+net.nₜ))]]

function ODE!(du, u, net, t)
	r = @view u[:,1]
	p = @view u[:,2]
	ϕ = @view u[:,3]
	ϕ .= max.(0., ϕ)
	du[:,1] .= drdt(net, r, p, ϕ)
	du[:,2] .= dpdt(net, r, p)
	du[1:net.nₚ+net.nₜ,3] .= dϕdt(net, view(p,1:net.nₚ+net.nₜ), view(ϕ,1:net.nₚ+net.nₜ))
end

steady_state_callback = TerminateSteadyState(1e-4, 1e-6)
default_duration = 24

"""
Save progression through time for plotting and inspection.
"""
function simulate(network::Network, duration::Integer=default_duration)
	problem = ODEProblem(ODE!, get_u₀(network), (0., duration*60.), network)
	solve(problem, callback=steady_state_callback, save_everystep=true)
end

"""
Find the steady-state solution.
"""
function steady_state(network::Network, u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	problem = ODEProblem(ODE!, u₀, (0., duration*60.), network)
	solve(problem, callback=steady_state_callback, save_everystep=false, save_start=false)
end
"""
Mutate a wildtype and get the steady state solution.
"""
function steady_state(network::Network, mutation::Union{<:AbstractVector,<:Int}, u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	steady_state(Network(network, mutation), u₀, duration)
end
function steady_states(network::Network, mutations::Matrix, u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	[steady_state(network, mutation, u₀, duration) for mutation in eachcol(mutations)]
end
function steady_states(network::Network, mutations, u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	[steady_state(network, mutation, u₀, duration) for mutation in mutations]
end

"""
The values we measure are log fold-change measurements of RNA expression levels in mutant compared to wildtype.
Note when we are using log or log2.
"""
logFC(wildtype, mutants::Vector) = hcat([log.(m[:,1,end] ./ wildtype[:,1,end]) for m in mutants]...)
"""
Default KOs is each TF and PK are knocked out in each individual experiment.
"""
function logFC(network::Network, mutations=1:network.nₜ+network.nₚ)
	u₀ = get_u₀(network)
	logFC(steady_state(network, u₀), steady_states(network, mutations, u₀))
end

end;
