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
@main function simulate(i=default_net, r="simulated_r.mat", p="simulated_p.mat", ϕ="simulated_phi.mat", t="simulated_t.mat")
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
@main function steadystate(mutate=nothing, i=default_net, r="steady_state_r" * (mutate == nothing ? ".mat" : "_$mutate.mat"), p="steady_state_p" * (mutate == nothing ? ".mat" : "_$mutate.mat"), ϕ="steady_state_phi" * (mutate == nothing ? ".mat" : "_$mutate.mat"); mutations=nothing)
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

@main function plot(i...; t="simulated_t.mat", o=stdout)
	simulations = [loaddlm(fname) for fname in i]
	time = loaddlm(t)
	p = Plotting.plot_simulation(time, simulations)
	Plotting.save(o, p)
	if o == stdout
		println("plotting complete")
		readline()
	end
	nothing
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

end;
