#!/usr/bin/env julia
"""
Defining and solving ODEs to the point of having the resulting simulated logFC values.
Used for simulation of gene expression levels.
"""
module ODEs
using DifferentialEquations: ODEProblem, solve, ODESolution, CallbackSet
using DiffEqCallbacks: TerminateSteadyState, PositiveDomain
using Distributions: Uniform
using LinearAlgebra: I
Main.@use "simulation/GeneRegulation"

export get_u₀
export timeseries, steady_state
export logFC
export @domainerror, timeseries_domain

const default_duration = 24 # hours
steady_state_callback = TerminateSteadyState(1e-4, 1e-6)
dtmax = 1  # minutes. used to limit step size to avoid instability. Reduction increases computational cost and stability.

"Helper function that uses a default value if given nothing, and otherwise uses what is given."
_dflt(::Nothing, default) = default
_dflt(value, ::Any) = value

""" ψ₀ == 0 for non-regulating proteins in order to have vectors match in length """
get_u₀(net::Network) = [net.r₀ net.p₀ [net.ψ₀; zeros(net.nₒ)]]
get_u₀(net::Network, r₀, p₀, ψ₀) = [_dflt(r₀,net.r₀) _dflt(p₀,net.p₀) [_dflt(ψ₀,net.ψ₀); zeros(net.nₒ)]]

function ODE!(du, u, net, t)
	u .= clamp.(u, 0., 1.) # can be used against dt <= dtmin warnings and domain errors.
    # Phos concentration cannot go above total concentration. 
    # It will naturally be enforced by dψdt except for certain values of λ₊, Wₚ₋, etc. 
    # so it is not guaranteed unless enforced here.
    u[:, 3] .= min.(u[:, 3], u[:, 2])
	r   = @view u[:, 1]
	p   = @view u[:, 2]
	ψₜₚ = @view u[1:net.nₜ+net.nₚ, 3]
	du[:, 1] .= drdt(net, r, p, ψₜₚ)
	du[:, 2] .= dpdt(net, r, p)
	du[1:net.nₜ+net.nₚ, 3] .= dψdt(net, view(p, 1:net.nₜ+net.nₚ), ψₜₚ)
end

default_callback(u₀) = CallbackSet(steady_state_callback, PositiveDomain(u₀))

"""
Save progression through time for plotting and inspection.
"""
function timeseries(network::Network; u₀::Matrix=get_u₀(network), duration::Real=default_duration)
	problem = ODEProblem(ODE!, u₀, (0., duration*60.), network)
	# it seems the callback is enforcing a save for every iteration so we get duplicate values with save_everystep=true.
	solve(problem, callback=default_callback(u₀), save_everystep=false, dtmax=dtmax, force_dtmin=true)
end
"""
Mutate a wildtype and simulate through time.
"""
timeseries(network::Network,  mutation, u₀, duration) = timeseries(Network(network, mutation); u₀=_dflt(u₀,get_u₀(network)), duration=_dflt(duration,default_duration))
timeseries(network::Network, ::Nothing, u₀, duration) = timeseries(network; u₀=_dflt(u₀,get_u₀(network)), duration=_dflt(duration,default_duration))
timeseries(network::Network,  mutation) = timeseries(Network(network, mutation))
timeseries(network::Network, ::Nothing) = timeseries(network)

"""
Find the steady-state solution.
"""
function steady_state(network::Network; u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	problem = ODEProblem(ODE!, u₀, (0., duration*60.), network)
	solve(problem, callback=default_callback(u₀), save_everystep=false, save_start=false, dtmax=dtmax, force_dtmin=true)
end
"""
Mutate a wildtype and get the steady state solution.
"""
function steady_state(network::Network, mutation::Union{<:AbstractVector,<:Integer}, u₀::Matrix=get_u₀(network), duration::Integer=default_duration)
	steady_state(Network(network, mutation); u₀=u₀, duration=duration)
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
	wildtype = steady_state(network; u₀=get_u₀(network))
	logFC(wildtype, steady_states(network, mutations, wildtype[:,:,end]))
end


"""
Call a function and catch any domain errors so we avoid the default verbose crash and instead display a non fatal error message and return nothing in that case.
return: macro returns the result of expression ex if there are no errors, nothing if there is a domain error or crashes if there is another error.
"""
macro domainerror(ex)
	quote try $(esc(ex)) catch e
		if e isa DomainError @error("Domain error: $(e.val)")
		else throw(e) end end
	end
end


timeseries_domain(args...) = @domainerror(timeseries(args...))
logFC_domain(args...) = @domainerror(logFC(args...))


end;
