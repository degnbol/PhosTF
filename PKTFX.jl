#!/usr/bin/env julia
include("simulation/GeneRegulation.jl")
include("simulation/ODEs.jl")
include("utilities/ReadWrite.jl")
include("utilities/CLI.jl")
include("Cytoscape.jl")
include("Plotting.jl")
include("ModelIteration.jl")
if !isdefined(Main, :Model) include("Model.jl") end # loaded by GeneRegulation
include("Weight.jl")
include("Inference.jl")
if !isdefined(Main, :ArrayUtils) include("utilities/ArrayUtils.jl") end


"Main caller for calling all project functions through command-line or julia REPL."
module PKTFX
using Fire
using Distributions: Uniform
using Plots
using ..ReadWrite, ..ArrayUtils
using ..Cytoscape, ..Plotting, ..ODEs, ..ModelIteration, ..Model, ..Weight
using ..GeneRegulation, ..Inference, ..CLI

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"


loadnet(i) = load(i, Network)

"""
Create random W.
"""
@main function W(n::Integer, nₜ::Integer=Integer(round(n/2 * .7)), nₚ::Integer=Integer(round(n/2 * .3)); WT_fname::String=default_Wₜ, WP_fname::String=default_Wₚ)
	savedlm(WT_fname, Weight.random_Wₜ(n , nₜ, nₚ))
	savedlm(WP_fname, Weight.random_Wₜ(nₜ, nₚ))
end

"""
- np: mandatory arg to set nₚ
"""
@main function correct(io=nothing, o=nothing; np=nothing)
	if np === nothing @error("supply --np"); return end
	i, o = inout(io, o)
	W = loaddlm(io)
	if Weight.correct!(W, np) @info("Corrections made.")
	else @info("NO corrections made.") end
	savedlm(o, W)
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

# load file(s) as a single 2D array regardless if they match in length along axis 1
hcatpad_load(fnames::Vector) = hcatpad(loaddlm(fname, Float64) for fname ∈ fnames)
hcatpad_load(fname::String) = loaddlm(fname, Float64)

"""
Write a graph defined by weight matrices to xgmml format.
"""
@main function xgmml(Wₜ, Wₚ; o=stdout, title=nothing, X=[])
	Wₜ, Wₚ = loaddlm(Wₜ), loaddlm(Wₚ)
	title = o == stdout ? "pktfx" : splitext(basename(o))[1]
	if isempty(X) write(o, Cytoscape.xgmml(Wₜ, Wₚ; title=title))
	else
		X = hcatpad_load(X)
		K = size(X,2)
		_,nₜ,nₚ = nₓnₜnₚ(Wₜ,Wₚ)
		# highlight each of the proteins if there are as many experiments as PKs+TFs
		highlight = nₜ+nₚ == K ? (1:K) : nothing
		write(o, Cytoscape.xgmml(Wₜ, Wₚ, X, highlight; title=title))
	end
end
"""
Write a graph defined by a simulation network file to xgmml format.
- X: node values. Each column of X is used for a separate copy of the graph.
"""
@main function xgmml(i=default_net; o=stdout, title=nothing, X=[])
	net = loadnet(i)
	title = o == stdout ? "pktfx" : splitext(basename(o))[1]
	if isempty(X) write(o, Cytoscape.xgmml(net; title=title))
	else
		X = hcatpad_load(X)
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
	Wₜ, Wₚ = loaddlm(Wₜ_fname), loaddlm(Wₚ_fname, Int64)
	n, nₜ = size(Wₜ)
	nₚ = size(Wₚ,2)
	@assert nₜ == size(Wₚ,1) - nₚ
	@assert n >= nₜ + nₚ
	save(o, Network(Wₜ, Wₚ))
end

@main function display(i=default_net; v::Integer=0)
	net = loadnet(i)
	println(net)
	if v > 0
		for g in net.genes println("\t", g)
			if v > 1
				for m in g.modules println("\t\t", m) end
			end
		end
	end
end

@main function estimateWt(i=default_net; o=stdout)
	Wₜ = GeneRegulation.estimate_Wₜ(loadnet(i))
	savedlm(o, Wₜ)
end


"""
Simulate a network.
- r,p,phi: fname. output files with time along axis 2 and protein along axis 1
- r0, p0, phi0: fname. starting values. Use this to let a mutated network start where the wildtype converged. 
Make them with either another simulate call, a steaady state call or write them manually.
"""
@main function simulate(mut_id=nothing, i=default_net; r=nothing, p=nothing, phi=nothing, t=nothing, duration=nothing, r0=nothing, p0=nothing, phi0=nothing)
	if   r  === nothing   r  = "sim_r"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   p  === nothing   p  = "sim_p"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if phi  === nothing phi  = "sim_phi" * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   t  === nothing   t  = "sim_t"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   r0 !== nothing   r0 = loaddlm(  r0)[:,end] end
	if   p0 !== nothing   p0 = loaddlm(  p0)[:,end] end
	if phi0 !== nothing phi0 = loaddlm(phi0)[:,end] end
	net = loadnet(i)
	u₀ = get_u₀(net, r0, p0, phi0)
	solution = @domainerror ODEs.simulate(net, mut_id, u₀, duration)
	if solution === nothing return end
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		savedlm(r,   solution[:,1,:])
		savedlm(p,   solution[:,2,:])
		savedlm(phi, solution[1:net.nₜ+net.nₚ,3,:])
		savedlm(t,   solution.t)
	end
end

"""
Plot simulations.
- nₚ: number of PK/PP.
- r,p,ϕ: fname. matrices with size=(#proteins, #times)
- t: fname. time vector
- o: optional file to write plot to
"""
@main function plot(nₚ::Integer, nₜ::Integer, r="sim_r.mat", p="sim_p.mat", ϕ="sim_phi.mat", t="sim_t.mat"; o=stdout)
	r, p, ϕ, t = loaddlm(r), loaddlm(p), loaddlm(ϕ), loaddlm(t)
	t = dropdims(t; dims=2)  # should be a column vector in file
	n = size(r,1); nₓ = n-(nₜ+nₚ)
	# protein along axis=1, mRNA,prot,phos along axis=2 and for measurements: time along axis=3
	values = zeros(n, 3, length(t))
	values[:,1,:], values[:,2,:], values[1:nₜ+nₚ,3,:] = r, p, ϕ
	nodes  = [["P$i" for i ∈ 1:nₚ]; ["T$i" for i ∈ 1:nₜ]; ["X$i" for i ∈ 1:nₓ]]
	labels = ["$i $l" for i ∈ nodes, l ∈ ["mRNA", "prot", "phos"]]
	styles = [s for i ∈ nodes, s ∈ [:solid, :solid, :dash]]
	widths = [w for i ∈ nodes, w ∈ [1, 2, 1]]
	# get [P, T, X] collections of data
	values = [reshape(values[1:nₚ,:,:], 3nₚ, :), reshape(values[nₚ+1:nₚ+nₜ,:,:], 3nₜ, :), values[nₚ+nₜ+1:end,1,:]]
	labels = [reshape(labels[1:nₚ,:], :), reshape(labels[nₚ+1:nₚ+nₜ,:], :), labels[nₚ+nₜ+1:end,1]]
	styles = [reshape(styles[1:nₚ,:], :), reshape(styles[nₚ+1:nₚ+nₜ,:], :), styles[nₚ+nₜ+1:end,1]]
	widths = [reshape(widths[1:nₚ,:], :), reshape(widths[nₚ+1:nₚ+nₜ,:], :), widths[nₚ+nₜ+1:end,1]]
	p = Plotting.plot_simulation(t, values, labels, styles, widths, ["PK/PP", "TF", "X"])
	Plotting.save(o, p)
	if o == stdout
		println("plotting complete")
		readline()
	end
	nothing
end


"""
Get the steady state levels for a network optionally mutating it.
- mut_id: index indicating which mutation to perform.
If "mutations" is not provided, it refers to index of the protein to mutate.
- i: network fname.
- r: out fname for r values.
- p: out fname for p values.
- ϕ: out fname for ϕ values.
- mut_file: optional fname to provided where "mut_id" will then be the index of a column.
If "mut" is not provided, the first (and ideally only) column of the file will be used.
"""
@main function steadystate(mut_id=nothing, i=default_net; r=nothing, p=nothing, phi=nothing, mut_file=nothing)
	if r   === nothing r   = "steady_r"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if p   === nothing p   = "steady_p"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if phi === nothing phi = "steady_phi" * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	net = loadnet(i)
	solution = steady_state(net, mut_id, mut_file)
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		savedlm(r,   solution[:,1,end])
		savedlm(p,   solution[:,2,end])
		savedlm(phi, solution[1:net.nₜ+net.nₚ,3,end])
	end
end
steady_state(net, ::Nothing, ::Nothing) = ODEs.steady_state(net)
steady_state(net, mutation::Integer, ::Nothing) = ODEs.steady_state(net, mutation)
steady_state(net, mutation_col, mutations_file::String) = steady_state(net, mutation_col, loaddlm(mutations_file, Int))
function steady_state(net, ::Nothing, mutations::Matrix)
	if size(mutations,2) == 1 return steady_state(net, 1, mutations)
	else error("column in mutation file not selected, use first arg") end
end
function steady_state(net, mutation_col::Integer, mutations::Matrix)
	if size(mutations,1) == net.nₚ+net.nₜ
		try mutations = convert(BitVector, mutations) catch InexactError end
	end
	ODEs.steady_state(net, mutations[:,mutation_col])
end

"""
Get the log fold-change values comparing mutant transcription levels to wildtype.
-	net: a network
-	o: stdout or file to write result to
"""
@main function logFC(net=default_net; o=stdout)
	measurements = @domainerror(ODEs.logFC(loadnet(net)))
	if measurements != nothing
		@info("logFC values simulated")
		savedlm(o, measurements)
	end
end
"""
Get the log fold-change values comparing mutant transcription levels to wildtype.
-	wt: wildtype expression level matrix. String fname
-	mut: mutant expression level matrix. String fname
-	muts...: additional mutant expression level matrices. list of string fname
-	o: stdout or file to write result to
"""
@main function logFC(wt, mut, muts...; o=stdout)
	wildtype = loaddlm(wt)
	mutants = [loaddlm(mutant) for mutant in [mut; muts...]]
	measurements = @domainerror ODEs.logFC(wildtype, mutants)
	if measurements != nothing
		@info("logFC values simulated")
		savedlm(o, measurements)
	end
end

"""
Get priors from files with the indicators 0=no edge, 1=possible edge, "+"=positive edge, "-"=negative edge.
Can be fed nothing values, and produces nothing values when a matrix would otherwise provide no additional information.
return: priors, priors_sign
"""
function _priors(WT_prior, WP_prior, n::Integer, nₜ::Integer, nₚ::Integer)
	if WT_prior === nothing && WP_prior === nothing return nothing, nothing end
	M, S = Model.priors(
	WT_prior === nothing ? n  : loaddlm(WT_prior),
	WP_prior === nothing ? nₚ : loaddlm(WP_prior))
	if all(Model._Wₜ(M,nₜ,nₚ) .== 1) && all(Model._Wₚ(M,nₜ,nₚ) .== 1) M = nothing end
	if all(S == 0) S = nothing end
	M, S
end

"""
Infer a weight matrix from logFC data.
- WT_prior/WP_prior: optionally limit Wₜ/Wₚ if they are partially known.
"""
@main function infer(X, nₜ::Integer, nₚ::Integer, ot="WT_infer.mat", op="WP_infer.mat"; WT_prior=nothing, WP_prior=nothing, epochs::Integer=5000, lambda::Real=1e-5)
	X = loaddlm(X, Float64)
	M, S = _priors(WT_prior, WP_prior, size(X,1), nₜ, nₚ)
	W = Inference.infer(X, nₜ, nₚ; epochs=epochs, λ=lambda, M=M, S=S)
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end
@main function inferB(B_LLC, nₚ::Integer, ot="WT_inferB.mat", op="WP_inferB.mat"; WT_prior=nothing, WP_prior=nothing, epochs::Integer=5000, lambda::Real=1e-5)
	B_LLC = loaddlm(B_LLC, Float64)
	nₚₜ = size(B_LLC,1); nₜ = nₚₜ-nₚ
	M, S = _priors(WT_prior, WP_prior, nₚₜ, nₜ, nₚ)
	W = Inference.infer_B(B_LLC, nₚ; epochs=epochs, λ=lambda, M=M, S=S)
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

"""
Remove edges less than a given threshold.
- thres: optional threshold
"""
@main function thres(io=nothing, o=nothing; thres=0.01)
	i, o = inout(io, o)
	mat = loaddlm(i, Float64)
	Weight.threshold!(mat, thres)
	savedlm(o, mat)
end


"""
Swap order of PK and TF in matrix (swap between PK-TF-X and TF-PK-X).
- n1: number of nodes in first group, which will be moved to become the second.
- n2: number of nodes in the second group, which will be move to become the first.
"""
@main function swapPT(io=nothing, o=nothing; n1=nothing, n2=nothing)
	i, o = inout(io, o)
	if n1 === nothing || n2 === nothing @error("Provide --n1 and --n2.") end
	mat = loaddlm(i, Float64)
	mat = reorder(mat, [n1+1:n1+n2;1:n1;n1+n2+1:maximum(size(mat))])
	savedlm(o, mat)
end

end;
