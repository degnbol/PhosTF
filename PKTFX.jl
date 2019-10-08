#!/usr/bin/env julia
include("GeneRegulation.jl")
"Main caller for calling all project functions through command-line or julia REPL."
module PKTFX
using Fire
using Distributions: Uniform
using Plots
include("utilities/ReadWrite.jl"); using .ReadWrite
include("utilities/CLI.jl")
include("Cytoscape.jl")
include("Plotting.jl")
include("ODEs.jl")
include("ModelIteration.jl")
include("Model.jl"); using .Model: nₓnₜnₚ
include("Weight.jl")

default_Wₜ = "WT.mat"
default_Wₚ = "WP.mat"
default_net = "net.bson"

"""
Create random W.
"""
@main function W(n::Integer, nₜ::Integer=Integer(round(n/2 * .7)), nₚ::Integer=Integer(round(n/2 * .3)); WT_fname::String=default_Wₜ, WP_fname::String=default_Wₚ)
	savedlm(WT_fname, Weight.Wₜ(n , nₜ, nₚ))
	savedlm(WP_fname, Weight.Wₜ(nₜ, nₚ))
end

@main function correct(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; ot="WT_cor.mat", op="WP_cor.mat")
	Wₜ, Wₚ = loaddlm(Wₜ_fname), loaddlm(Wₚ_fname)
	if Weight.correct!(Wₜ, Wₚ) @info("Corrections made.") end
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

@main function iteratemodel(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

@main function xgmml(Wₜ="WT.mat", Wₚ="WP.mat"; X=nothing, o=stdout)
	Wₜ, Wₚ = loaddlm(Wₜ), loaddlm(Wₚ)
	if X != nothing X = loaddlm(X) end
	_,nₜ,nₚ = nₓnₜnₚ(Wₜ,Wₚ)
	K = size(X,2)
	# highlight each of the proteins if there are as many experiments as PKs+TFs
	highlight = nₜ+nₚ == K ? (1:K) : nothing
	write(o, Cytoscape.xgmml(Wₜ, Wₚ, X, highlight))
end


"""
Create a random network from W.
"""
@main function network(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; o::String=default_net)
	save(o, Main.GeneRegulation.Network(loaddlm(Wₜ_fname), loaddlm(Wₚ_fname)))
end

@main function display(i=default_net; v::Integer=0)
	net = load(i, Main.GeneRegulation.Network)
	println(net)
	if v > 0 for g in net.genes
		println("\t", g)
		if v > 1 for m in g.modules
			println("\t\t", m)
		end end
	end end
end

@main function estimateWt(i=default_net; o=stdout)
	net = load(i, Main.GeneRegulation.Network)
	Wₜ = Main.GeneRegulation.estimate_Wₜ(net)
	savedlm(o, Wₜ)
end


"""
Simulate a network.
"""
@main function simulate(i=default_net, r="sim_r.mat", p="sim_p.mat", ϕ="sim_phi.mat", t="sim_t.mat")
	net = load(i, Main.GeneRegulation.Network)
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
	net = load(i, Main.GeneRegulation.Network)
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
		net = load(i[1], Main.GeneRegulation.Network)
		measurements = ODEs.logFC(net)
	else
		wildtype = loaddlm(i[1])
		mutants = [loaddlm(mutant) for mutant in i[2:end]]
		measurements = ODEs.logFC(wildtype, mutants)
	end
	savedlm(o, measurements)
end


end;
