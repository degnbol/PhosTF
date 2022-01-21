#!/usr/bin/env julia
SRC = readchomp(`git root`) * "/src/"
isdefined(Main, :GeneRegulation) || include("GeneRegulation.jl")
isdefined(Main, :ODEs) || include(SRC * "simulation/ODEs.jl")
isdefined(Main, :ReadWrite) || include(SRC * "utilities/ReadWrite.jl")
using Fire
using LinearAlgebra
using Plots
using Test  # @test_logs used in randomNetLogFCs
using .Threads: @threads
using .ReadWrite

# defaults
default_Wₜ, default_Wₚ = "WT.adj", "WP.adj"
default_net = "net.bson"


loadnet(i) = ReadWrite.load(i, GeneRegulation.Network)

"""
Create a random network from Wₜ and Wₚ.
- Wₜfname, Wₚfname: filenames of Wₜ and Wₚ matrices with edge presence and sign integer data. 
    Sizes are nᵥ×nₜ and nₜ+nₚ×nₚ, and the first dimension should be sorted TF, KP, O.
- header: bool indicating if weight matrices contains headers (and then possibly row names).
    Note that it is not (yet) implemented to store this info in the generated network model.
"""
@main function network(Wₜfname::String=default_Wₜ, Wₚfname::String=default_Wₚ; header::Bool=false, o::String=default_net)
	Wₜ = loaddlm(Wₜfname; header=header)
	# make sure to enforce that it is Int which indicates weight presence and sign as opposed to Float that indicates weight magnitude.
	Wₚ = loaddlm(Wₚfname, Int; header=header)
    if header
        # then there should also be a column with rownames
        @assert eltype(Wₜ[!, 1]) <: AbstractString "Wₜ has a header but no row names."
        @assert eltype(Wₚ[!, 1]) <: AbstractString "Wₚ has a header but no row names."
        TFs = names(Wₜ)[2:end]
        KPs = names(Wₚ)[2:end]
        TFKPs = [TFs; KPs]
        # assume order is the same in column and row names. Could be changed to reorder.
        @assert all(Wₜ[1:length(TFKPs), 1] .== TFKPs)
        @assert all(Wₚ[!, 1] .== TFKPs)
        Wₜ = Matrix(Wₜ[:, 2:end])
        Wₚ = Matrix(Wₚ[:, 2:end])
    end
	nᵥ, nₜ = size(Wₜ)
	nₚ = size(Wₚ, 2)
    @assert nₜ+nₚ == size(Wₚ, 1) "nₜ+nₚ != size(Wₚ, 1). nₜ=$nₜ. nₚ=$nₚ."
	@assert nᵥ >= nₜ + nₚ
	save(o, GeneRegulation.Network(Wₜ, Wₚ))
end

@main function printNet(i=default_net; v::Integer=0)
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
Simulate a network and return the full time series.
- r,p,psi: fname. output files with time along axis 2 and protein along axis 1
- r0, p0, psi0: fname. starting values. Use this to let a mutated network start where the wildtype converged. 
Make them with either another simulate call, a steaady state call or write them manually.
"""
@main function timeseries(mut_id=nothing, i=default_net; r=nothing, p=nothing, psi=nothing, t=nothing, duration=nothing, r0=nothing, p0=nothing, psi0=nothing)
    # default output names
    mutSuf = mut_id === nothing ? "" : "_$mut_id"
	if   r  === nothing   r  = "sim_r$mutSuf.mat" end
	if   p  === nothing   p  = "sim_p$mutSuf.mat" end
	if psi  === nothing psi  = "sim_psi$mutSuf.mat" end
	if   t  === nothing   t  = "sim_t$mutSuf.mat" end
    
	if   r0 !== nothing   r0 = ReadWrite.loaddlm(  r0)[end, :] end
	if   p0 !== nothing   p0 = ReadWrite.loaddlm(  p0)[end, :] end
	if psi0 !== nothing psi0 = ReadWrite.loaddlm(psi0)[end, :] end
	
    net = loadnet(i)
	u₀ = ODEs.get_u₀(net, r0, p0, psi0)
    @assert !any(isnan.(u₀))
	solution = ODEs.@domainerror ODEs.timeseries(net, mut_id, u₀, duration)
	solution === nothing && return
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
        # dimensions are (node, RNA/protein/activation, time)
        # transpose to have time along first dim.
		savedlm(r,   solution[:, 1, :]')
		savedlm(p,   solution[:, 2, :]')
		savedlm(psi, solution[1:net.nₜ+net.nₚ, 3, :]')
        savedlm(t,   solution.t)
	end
end

"""
Plot time series.
- nₜ: number of TFs.
- nₚ: number of KPs.
- r,p,ψ: fname. matrices with size=(#proteins, #times)
- t: fname. time vector
- o: optional file to write plot to
"""
@main function plot(nₜ::Integer, nₚ::Integer, r="sim_r.mat", p="sim_p.mat", ψ="sim_psi.mat", t="sim_t.mat"; o=stdout)
    include(SRC * "simulation/PlotTimeSeries.jl") # save startup time by only including it when necessary

	# shape is times × nodes for r, p, ψ
    r, p, ψ, t = ReadWrite.loaddlm.([r, p, ψ, t])
	t = dropdims(t; dims=2)  # should be a column vector in file
    # dimensions should be time points along first dim, nodes along second in all files.
	nᵥ = size(r, 2)
    nₒ = nᵥ - (nₜ + nₚ)

	values = zeros(nᵥ, 3, length(t))
	values[:, 1, :], values[:, 2, :], values[1:nₜ+nₚ, 3, :] = r', p', ψ'
    meass = ["mRNA", "prot", "phos"]
	nodes  = [["KP$i" for i ∈ 1:nₚ]; ["TF$i" for i ∈ 1:nₜ]; ["O$i" for i ∈ 1:nₒ]]
	labels = ["$i $l" for i ∈ nodes, l ∈ meass]
	styles = [s for i ∈ nodes, s ∈ [:solid, :solid, :dash]]
	widths = [w for i ∈ nodes, w ∈ [1, 2, 1]]
    colors = [c for c ∈ ["#B663B1", "#6DB845", "#545000"], l ∈ meass]
	# get [KP, TF, O] collections of data
	values = [reshape(values[1:nₚ,:,:], 3nₚ, :), reshape(values[nₚ+1:nₚ+nₜ,:,:], 3nₜ, :), values[nₚ+nₜ+1:end,1,:]]
	labels = [reshape(labels[1:nₚ,:], :), reshape(labels[nₚ+1:nₚ+nₜ,:], :), labels[nₚ+nₜ+1:end,1]]
	styles = [reshape(styles[1:nₚ,:], :), reshape(styles[nₚ+1:nₚ+nₜ,:], :), styles[nₚ+nₜ+1:end,1]]
	widths = [reshape(widths[1:nₚ,:], :), reshape(widths[nₚ+1:nₚ+nₜ,:], :), widths[nₚ+nₜ+1:end,1]]
	p = PlotTimeSeries.plot_timeseries(t, values, labels, styles, widths, ["KP", "TF", "O"]; colors=colors)
	if o != stdout
        savefig(p, o)
    else
		println("plotting complete")
        readline() # enter will end program
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
		savedlm(r,   solution[:, 1, end])
		savedlm(p,   solution[:, 2, end])
		savedlm(psi, solution[1:net.nₜ+net.nₚ, 3, end])
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
- net: a network
- o: stdout or file to write result to
"""
@main function logFC(net=default_net; o=stdout)
    net = loadnet(net)
	measurements = ODEs.@domainerror(ODEs.logFC(net))
	if measurements !== nothing
		@info("logFC values simulated")
        # add some default names to the logFC table, custom names not implemented (yet).
        names = vcat(["TF$i" for i in 1:net.nₜ], ["KP$i" for i in 1:net.nₚ], ["O$i" for i in 1:net.nₒ])
        savedlm(o, measurements; colnames=names[1:net.nₜ+net.nₚ], rownames=names)
	end
end
"""
Get the log fold-change values comparing mutant transcription levels to wildtype.
- wt: wildtype expression level matrix. String fname
- mut: mutant expression level matrix. String fname
- muts...: additional mutant expression level matrices. list of string fname
- o: stdout or file to write result to
"""
@main function logFC(wt, mut, muts...; o=stdout)
	wildtype = ReadWrite.loaddlm(wt)
	mutants = [ReadWrite.loaddlm(mutant) for mutant in [mut; muts...]]
	measurements = ODEs.@domainerror(ODEs.logFC(wildtype, mutants))
	if measurements !== nothing
		@info("logFC values simulated")
		savedlm(o, measurements)
	end
end

"""
Given an adjacency matrix without KP, TF, O assignment, generate a random network and simulate it, returning the logFC values.
- B: path of gold standard matrix which is an adjacency matrix without assignments of node types.
- dir: path of output directory to store result files. All files are saved with default names, e.g. WT.mat, WP.mat, net.bson, X_sim.mat
- max_np: Try to get this many nodes assigned as KP. Upper limit. Default=30% of nodes.
- max_attempts: Maximum number of tries to generate a random net before giving up if there keeps being a convergence or other error from simulation.
"""
@main function randomNetLogFCs(B::String, dir::String; max_np::Union{Int,Nothing}=nothing, max_attempts::Int=3)
    # suppress Fire output with ;
    include(SRC * "/src/simulation/weight.jl");  # random(...) and correct()
    # set max_nₚ to default 30% of nodes if not specified.
    if max_np === nothing
        nᵥ = size(ReadWrite.loaddlm(B), 1)
        max_np = ceil(Int, nᵥ * .3)
    end

    println(dir)
    Wₜfname = joinpath(dir, default_Wₜ)
    Wₚfname = joinpath(dir, default_Wₚ)
    netfname = joinpath(dir, default_net)
    logFCfname = joinpath(dir, "X_sim.mat")
    
    for attempt in 1:max_attempts
        attempt == 1 || println("$dir attempt $attempt")
        # make WT and WP with {-1, 0, 1} for deactvation, no edge and activation
        random(B, max_np; WT=Wₜfname, WP=Wₚfname)
        # there should be nothing to correct
        #= correct(Wₜfname, Wₚfname; ot=Wₜfname, op=Wₚfname) =#
        # make net.bson which has a fully defined network instance with all it's constants, etc.
        network(Wₜfname, Wₚfname; o=netfname)

        rm(logFCfname; force=true)
        # simulate logFC values and test that there are no warnings
        try
            @test_logs (:info, "$dir logFC values simulated") logFC(netfname; o=logFCfname)
        catch
            continue
        end
        break
    end
end

