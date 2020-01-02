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
using LinearAlgebra
using Plots
using ..ReadWrite, ..ArrayUtils
using ..Cytoscape, ..Plotting, ..ODEs, ..ModelIteration, ..Model, ..Weight
using ..GeneRegulation, ..Inference, ..CLI

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"


loadnet(i) = load(i, Network)

"""
Create random W from a adjacency matrix containing B.
"""
@main function W(B::String, nₚₖ::Integer, nₚₚ::Integer; WT_fname::String=default_Wₜ, WP_fname::String=default_Wₚ)
	Wₜ, Wₚ = Weight.random_W(loaddlm(B), nₚₖ, nₚₚ)
	savedlm(WT_fname, Wₜ)
	savedlm(WP_fname, Wₚ)
end

"""
- np: mandatory arg to set nₚ
- save: should we save files if no corrections are made?
"""
@main function correct(io=nothing, o=nothing; np=nothing, save::Bool=false)
	if np === nothing @error("supply --np"); return end
	i, o = inout(io, o)
	W = loaddlm(io)
	if Weight.correct!(W, np)
		@info("Corrections made.")
		savedlm(o, W)
	else
		@info("NO corrections made.")
		if save savedlm(o, W) end
	end
end
@main function correct(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; ot="WT_cor.mat", op="WP_cor.mat", save::Bool=false)
	Wₜ, Wₚ = loaddlm(Wₜ_fname), loaddlm(Wₚ_fname)
	if Weight.correct!(Wₜ, Wₚ)
		@info("Corrections made.")
		savedlm(ot, Wₜ)
		savedlm(op, Wₚ)
	else
		@info("NO corrections made.")
		if save
			savedlm(ot, Wₜ)
			savedlm(op, Wₚ)
		end
	end
end

@main function iteratemodel(Wₜ="WT.mat", Wₚ="WP.mat"; o=stdout)
	X = ModelIteration.converge(loaddlm(Wₜ), loaddlm(Wₚ))
	savedlm(o, X)
end

"Load file(s) as a single 2D array regardless if they match in length along axis 1."
hcatpad_load(fnames::Vector) = hcatpad(loaddlm(fname, Float64) for fname ∈ fnames)
hcatpad_load(fname::String) = loaddlm(fname, Float64)

"""
Write a graph defined by weight matrices to xgmml format.
- X: node values. Each column of X is used for a separate copy of the graph.
"""
@main function xgmml(Wₜ, Wₚ::String; o=stdout, title=nothing, X=[])
	Wₜ, Wₚ = loaddlm(Wₜ), loaddlm(Wₚ)
	_,nₜ,nₚ = nₓnₜnₚ(Wₜ,Wₚ)
	xgmml([Wₜ, Wₚ], o, nₜ, nₚ, title, X)
end
@main function xgmml(i=default_net; o=stdout, title=nothing, X=[])
	net = loadnet(i)
	xgmml([net], o, net.nₜ, net.nₚ, title, X)
end
@main function xgmml(W, nₜ::Integer, nₚ::Integer; o=stdout, title=nothing, X=[])
	W = loaddlm(W)
	Wₜ, Wₚ = W[:,nₚ+1:nₚ+nₜ], W[:,1:nₚ] # we don't use the Model.WₜWₚ function since we want to allow P→X edges in case W==T
	xgmml([Wₜ, Wₚ], o, nₜ, nₚ, title, X)
end
function xgmml(i, o, nₜ, nₚ, title=nothing, X=[])
	if title === nothing title = o == stdout ? "pktfx" : splitext(basename(o))[1] end
	if isempty(X) write(o, Cytoscape.xgmml(i...; title=title))
	else
		X = hcatpad_load(X)
		K = size(X,2)
		# highlight each of the proteins if there are as many experiments as PKs+TFs
		highlight = nₜ+nₚ == K ? (1:K) : nothing
		write(o, Cytoscape.xgmml(i..., X, highlight; title=title))
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

"Get the Wₚ matrix that contains edges indicating activation as opposed to phosphorylation since we don't infer phosphorylation but activation."
@main function WP_activation(i=default_net; o=stdout)
	net = loadnet(i)
	Wₚ_activation = GeneRegulation.Wₚ_activation(net)
	Wₚ = sign.(net.Wₚₖ .- net.Wₚₚ)
	phos_regulated = sum(any(Wₚ .!= 0; dims=2))
	changes = sum(Wₚ .!= Wₚ_activation)
	@info("""$changes edges were changed.
	$(sum(net.phos_activation)) out of $(length(net.phos_activation)) TFs/PKs/PPs are phos activated.
	$phos_regulated are phos regulated.""")
	savedlm(o, Wₚ_activation)
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
- lambdaWT: whether the WT edges are regularized. Should only be used if highly trusted WT_priors are provided.
- PKPP: fname. vector, each element is -1 or 1 indicating PP or PK. 0s ignored.
- WT/WP: previous run to continue.
"""
@main function infer(X, nₜ::Integer, nₚ::Integer, ot="WT_infer.mat", op="WP_infer.mat"; epochs::Integer=5000, 
	lambda::Real=.1, lambdaW::Real=0., lambdaWT::Bool=true, WT_prior=nothing, WP_prior=nothing, PKPP=nothing, WT=nothing, WP=nothing, J=nothing)
	ot, op = abspath(ot), abspath(op)  # workaround for weird cwd issues
	X = loaddlm(X, Float64)
	if J !== nothing
		J = loaddlm(J, Float64)
		@assert size(J) == size(X)
	end
	n = size(X,1)
	M, S = _priors(WT_prior, WP_prior, n, nₜ, nₚ)
	
	if PKPP !== nothing
		PKPP = vec(loaddlm(PKPP))
		padding = zeros(n-length(PKPP))
		Iₚₖ = diagm([PKPP .== +1; padding])
		Iₚₚ = diagm([PKPP .== -1; padding])
	else
		Iₚₖ, Iₚₚ = nothing, nothing
	end

	W = (WT === nothing || WP === nothing) ? nothing : Model._W(loaddlm(WT), loaddlm(WP))

	W = Inference.infer(X, nₜ, nₚ; epochs=epochs, λ=lambda, λW=lambdaW, λWT=lambdaWT, M=M, S=S, Iₚₖ=Iₚₖ, Iₚₚ=Iₚₚ, W=W, J=J)
	Wₜ, Wₚ = Model.WₜWₚ(W, nₜ, nₚ)
	savedlm(ot, Wₜ)
	savedlm(op, Wₚ)
end

"""
Remove edges less than a given threshold.
- thres: optional threshold
"""
@main function thres(io=nothing, o=nothing; thres=0.001)
	i, o = inout(io, o)
	mat = loaddlm(i, Float64)
	Weight.threshold!(mat, thres)
	savedlm(o, mat)
end

"""
Get the total effects T from single KO experiments.
"""
@main function T(io=nothing, o=nothing)
	i, o = inout(io, o)
	X = loaddlm(i)
	savedlm(o, Model.X2T(X))
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

"""
input: WP.mat
output: a vector indicating PK with -1, and PP with 1. 
"""
@main function PKPP(io=nothing, o=nothing)
	i, o = inout(io, o)
	savedlm(o, Weight.PKPP(loaddlm(i)))
end

end;
