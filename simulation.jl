#!/usr/bin/env julia
include("src/simulation/GeneRegulation.jl")
isdefined(Main, :ODEs) || include("src/simulation/ODEs.jl")
include("src/utilities/ReadWrite.jl")
(length(ARGS) == 0 || ARGS[1] == "plot") && include("src/Plotting.jl")


using Fire
using LinearAlgebra
using Plots

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"


loadnet(i) = ReadWrite.load(i, GeneRegulation.Network)

"""
Create a random network from W.
"""
@main function network(Wₜ_fname::String=default_Wₜ, Wₚ_fname::String=default_Wₚ; o::String=default_net)
	Wₜ, Wₚ = ReadWrite.loaddlm(Wₜ_fname), ReadWrite.loaddlm(Wₚ_fname, Int)
	nᵥ, nₜ = size(Wₜ)
	nₚ = size(Wₚ,2)
	@assert nₜ == size(Wₚ,1) - nₚ
	@assert nᵥ >= nₜ + nₚ
	save(o, GeneRegulation.Network(Wₜ, Wₚ))
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
- r,p,psi: fname. output files with time along axis 2 and protein along axis 1
- r0, p0, psi0: fname. starting values. Use this to let a mutated network start where the wildtype converged. 
Make them with either another simulate call, a steaady state call or write them manually.
"""
@main function simulate(mut_id=nothing, i=default_net; r=nothing, p=nothing, psi=nothing, t=nothing, duration=nothing, r0=nothing, p0=nothing, psi0=nothing)
	if   r  === nothing   r  = "sim_r"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   p  === nothing   p  = "sim_p"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if psi  === nothing psi  = "sim_psi" * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   t  === nothing   t  = "sim_t"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if   r0 !== nothing   r0 = ReadWrite.loaddlm(  r0)[:,end] end
	if   p0 !== nothing   p0 = ReadWrite.loaddlm(  p0)[:,end] end
	if psi0 !== nothing psi0 = ReadWrite.loaddlm(psi0)[:,end] end
	net = loadnet(i)
	u₀ = ODEs.get_u₀(net, r0, p0, psi0)
    @assert !any(isnan.(u₀))
	solution = ODEs.@domainerror ODEs.simulate(net, mut_id, u₀, duration)
	solution === nothing && return
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		ReadWrite.savedlm(r,   solution[:,1,:])
		ReadWrite.savedlm(p,   solution[:,2,:])
		ReadWrite.savedlm(psi, solution[1:net.nₜ+net.nₚ,3,:])
		ReadWrite.savedlm(t,   solution.t)
	end
end

"""
Plot simulations.
- nₚ: number of KPs.
- r,p,ψ: fname. matrices with size=(#proteins, #times)
- t: fname. time vector
- o: optional file to write plot to
"""
@main function plot(nₚ::Integer, nₜ::Integer, r="sim_r.mat", p="sim_p.mat", ψ="sim_psi.mat", t="sim_t.mat"; o=stdout)
    r, p, ψ, t = ReadWrite.loaddlm.([r, p, ψ, t])
	t = dropdims(t; dims=2)  # should be a column vector in file
	nᵥ = size(r,1); nₒ = nᵥ-(nₜ+nₚ)
	# protein along axis=1, mRNA,prot,phos along axis=2 and for measurements: time along axis=3
	values = zeros(nᵥ, 3, length(t))
	values[:,1,:], values[:,2,:], values[1:nₜ+nₚ,3,:] = r, p, ψ
	nodes  = [["P$i" for i ∈ 1:nₚ]; ["T$i" for i ∈ 1:nₜ]; ["O$i" for i ∈ 1:nₒ]]
	labels = ["$i $l" for i ∈ nodes, l ∈ ["mRNA", "prot", "phos"]]
	styles = [s for i ∈ nodes, s ∈ [:solid, :solid, :dash]]
	widths = [w for i ∈ nodes, w ∈ [1, 2, 1]]
	# get [KP, TF, O] collections of data
	values = [reshape(values[1:nₚ,:,:], 3nₚ, :), reshape(values[nₚ+1:nₚ+nₜ,:,:], 3nₜ, :), values[nₚ+nₜ+1:end,1,:]]
	labels = [reshape(labels[1:nₚ,:], :), reshape(labels[nₚ+1:nₚ+nₜ,:], :), labels[nₚ+nₜ+1:end,1]]
	styles = [reshape(styles[1:nₚ,:], :), reshape(styles[nₚ+1:nₚ+nₜ,:], :), styles[nₚ+nₜ+1:end,1]]
	widths = [reshape(widths[1:nₚ,:], :), reshape(widths[nₚ+1:nₚ+nₜ,:], :), widths[nₚ+nₜ+1:end,1]]
	p = Plotting.plot_simulation(t, values, labels, styles, widths, ["KP", "TF", "O"])
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
- psi: out fname for ψ values.
- mut_file: optional fname to provided where "mut_id" will then be the index of a column.
If "mut" is not provided, the first (and ideally only) column of the file will be used.
"""
@main function steadystate(mut_id=nothing, i=default_net; r=nothing, p=nothing, psi=nothing, mut_file=nothing)
	if r   === nothing r   = "steady_r"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if p   === nothing p   = "steady_p"   * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	if psi === nothing psi = "steady_psi" * (mut_id === nothing ? "" : "_$mut_id") * ".mat" end
	net = loadnet(i)
	solution = steady_state(net, mut_id, mut_file)
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		ReadWrite.savedlm(r,   solution[:,1,end])
		ReadWrite.savedlm(p,   solution[:,2,end])
		ReadWrite.savedlm(psi, solution[1:net.nₜ+net.nₚ,3,end])
	end
end
steady_state(net, ::Nothing, ::Nothing) = ODEs.steady_state(net)
steady_state(net, mutation::Integer, ::Nothing) = ODEs.steady_state(net, mutation)
steady_state(net, mutation_col, mutations_file::String) = steady_state(net, mutation_col, ReadWrite.loaddlm(mutations_file, Int))
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
	measurements = ODEs.@domainerror(ODEs.logFC(loadnet(net)))
	if measurements !== nothing
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
	wildtype = ReadWrite.loaddlm(wt)
	mutants = [ReadWrite.loaddlm(mutant) for mutant in [mut; muts...]]
	measurements = ODEs.@domainerror(ODEs.logFC(wildtype, mutants))
	if measurements !== nothing
		@info("logFC values simulated")
		ReadWrite.savedlm(o, measurements)
	end
end

