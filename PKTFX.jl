#!/usr/bin/env julia
include("simulation/GeneRegulation.jl")
include("simulation/ODEs.jl")
include("utilities/ReadWrite.jl")
include("utilities/CLI.jl")
include("Cytoscape.jl")
include("Plotting.jl")
include("ModelIteration.jl")
include("Model.jl")
include("Weight.jl")
include("Inference.jl")


"Main caller for calling all project functions through command-line or julia REPL."
module PKTFX
using Fire
using Distributions: Uniform
using Plots
using ..ReadWrite
using ..Cytoscape, ..Plotting, ..ODEs, ..ModelIteration, ..Model, ..Weight
using ..GeneRegulation, ..Inference, ..CLI


default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"

loadnet(i) = load(i, Network)

"""
Create random W.
"""
@main function W(n::Integer, nₜ::Integer=Integer(round(n/2 * .7)), nₚ::Integer=Integer(round(n/2 * .3)); WT_fname::String=default_Wₜ, WP_fname::String=default_Wₚ)
	savedlm(WT_fname, Weight.Wₜ(n , nₜ, nₚ))
	savedlm(WP_fname, Weight.Wₜ(nₜ, nₚ))
end

@main function correct(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; ot="WT_cor.mat", op="WP_cor.mat")
	Wₜ, Wₚ = loaddlm(Wₜ_fname), loaddlm(Wₚ_fname)
	if Weight.correct!(Wₜ, Wₚ) @info("Corrections made.")
	else @info("NO corrections made.") end
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

@main function iteratemodel(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

"""
Write a graph defined by weight matrices to xgmml format.
"""
@main function xgmml(Wₜ, Wₚ; o=stdout, X=nothing, title=nothing)
	Wₜ, Wₚ = loaddlm(Wₜ), loaddlm(Wₚ)
	title = o == stdout ? "pktfx" : splitext(basename(o))[1]
	if X == nothing write(o, Cytoscape.xgmml(Wₜ, Wₚ; title=title))
	else
		X = loaddlm(X, Float64)
		K = size(X,2)
		_,nₜ,nₚ = nₓnₜnₚ(Wₜ,Wₚ)
		# highlight each of the proteins if there are as many experiments as PKs+TFs
		highlight = nₜ+nₚ == K ? (1:K) : nothing
		write(o, Cytoscape.xgmml(Wₜ, Wₚ, X, highlight; title=title))
	end
end
"""
Write a graph defined by a simulation network file to xgmml format.
"""
@main function xgmml(i=default_net; o=stdout, X=nothing, title=nothing)
	net = loadnet(i)
	title = o == stdout ? "pktfx" : splitext(basename(o))[1]
	if X == nothing write(o, Cytoscape.xgmml(net; title=title))
	else
		X = loaddlm(X, Float64)
		K = size(X,2)
		# highlight each of the proteins if there are as many experiments as PKs+TFs
		highlight = net.nₜ+net.nₚ == K ? (1:K) : nothing
		write(o, Cytoscape.xgmml(net, X, highlight; title=title))
	end
end
"""
Create a random network from W.
"""
@main function network(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; o::String=default_net)
	save(o, Network(loaddlm(Wₜ_fname), loaddlm(Wₚ_fname)))
end

@main function display(i=default_net; v::Integer=0)
	net = loadnet(i)
	println(net)
	if v > 0 for g in net.genes
		println("\t", g)
		if v > 1 for m in g.modules
			println("\t\t", m)
		end end
	end end
end

@main function estimateWt(i=default_net; o=stdout)
	Wₜ = GeneRegulation.estimate_Wₜ(loadnet(i))
	savedlm(o, Wₜ)
end


"""
Simulate a network.
"""
@main function simulate(i=default_net, r="sim_r.mat", p="sim_p.mat", ϕ="sim_phi.mat", t="sim_t.mat")
	net = loadnet(i)
	solution = ODEs.simulate(net)
	println(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		savedlm(r, solution[:,1,:])
		savedlm(p, solution[:,2,:])
		savedlm(ϕ, solution[1:net.nₜ+net.nₚ,3,:])
		savedlm(t, solution.t)
	end
end

"""
Plot simulations.
- i: matrices with size=(#proteins, #times)
- o: optional file to write plot to
- ids: list of protein ids to only plot those.
"""
@main function plot(i...; t="sim_t.mat", o=stdout, ids=[])
	simulations = [loaddlm(fname) for fname in i]
	if !isempty(ids)
		simulations = [sim[ids,:] for sim in simulations]
	end
	time = loaddlm(t)
	p = Plotting.plot_simulation(time, simulations)
	Plotting.save(o, p)
	if o == stdout
		println("plotting complete")
		readline()
	end
	nothing
end


"""
Get the steady state levels for a network optionally mutating it.
- mutate: index indicating which mutation to perform.
If "mutations" is not provided, it refers to index of the protein to mutate.
- i: network fname.
- r: out fname for r values.
- p: out fname for p values.
- ϕ: out fname for ϕ values.
- mutations: optional fname to provided where "mutate" will then be the index of a column.
If "mutate" is not provided, the first (and ideally only) column of the file will be used.
"""
@main function steadystate(mutate=nothing, i=default_net, r="steady_r" * (mutate == nothing ? ".mat" : "_$mut.mat"), p="steady_p" * (mutate == nothing ? ".mat" : "_$mut.mat"), ϕ="steady_phi" * (mutate == nothing ? ".mat" : "_$mut.mat"); mutations=nothing)
	net = loadnet(i)
	if mutations != nothing
		if mutate == nothing mutate = 1 end
		mutations = loaddlm(mutations, Int)[:,mutate]
		if length(mutations) == net.nₚ+net.nₜ
			try mutations = convert(BitVector, mutations) catch InexactError end
		end
		solution = ODEs.steady_state(net, mutations)
	else
		solution = mutate == nothing ? ODEs.steady_state(net) : ODEs.steady_state(net, mutate)
	end
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		savedlm(r, solution[:,1,end])
		savedlm(p, solution[:,2,end])
		savedlm(ϕ, solution[1:net.nₜ+net.nₚ,3,end])
	end
end

"""
Get the log fold-change values comparing mutant transcription levels to wildtype.
-	i: either a network or a list of expression levels where the first one listed is the wildtype
-	o: stdout or file to write result to
"""
@main function logFC(i...; o=stdout)
	if length(i) == 1
		measurements = ODEs.logFC(loadnet(i[1]))
	else
		wildtype = loaddlm(i[1])
		mutants = [loaddlm(mutant) for mutant in i[2:end]]
		measurements = ODEs.logFC(wildtype, mutants)
	end
	@info("logFC values simulated")
	savedlm(o, measurements)
end

"""
Infer a weight matrix from logFC data.
"""
@main function infer(X, nₜ::Integer, nₚ::Integer, ot="WT_infer.mat", op="WP_infer.mat"; epochs::Integer=10000, lambda::AbstractFloat=.1)
	W = Inference.infer(loaddlm(X, Float64), nₜ, nₚ; epochs=epochs, λ=lambda)
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

@main function thres(io=nothing, o=nothing; thres=0.01)
	i, o = inout(io, o)
	mat = loaddlm(i)
	Weight.threshold!(mat, thres)
	savedlm(o, mat)
end

end;
