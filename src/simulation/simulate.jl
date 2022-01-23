#!/usr/bin/env julia
@src "simulation/GeneRegulation"
@src "simulation/ODEs"
@use "utilities/ReadWrite"
using Fire
using LinearAlgebra
using Plots
using Test  # @test_logs used in randomNetLogFCs
using .Threads: @threads

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
Read timeseries from file and optionally assert the names match those in a network.
"""
function read_timeseries(i::String, net=nothing)
    df = ReadWrite.loaddlm(i; header=true)
    mat, colnames = Matrix(df), names(df)
    r0idx = startswith.(colnames, 'r')
    p0idx = startswith.(colnames, 'p')
    ψ0idx = startswith.(colnames, 'ψ')
    t = df[:, "t"] 
    r = mat[:, r0idx]
    p = mat[:, p0idx]
    ψ = mat[:, ψ0idx]
    if net !== nothing
        @assert all(chop.(colnames[r0idx], head=2, tail=0) .== net.names) "$([s[3:end] for s in colnames[r0idx]]) != $(net.names)"
        @assert all(chop.(colnames[p0idx], head=2, tail=0) .== net.names) "$([s[3:end] for s in colnames[p0idx]]) != $(net.names)"
        @assert all(chop.(colnames[ψ0idx], head=2, tail=0) .== net.names[1:net.nₜ+net.nₚ]) "$([s[3:end] for s in colnames[ψ0idx]]) != $(net.names[1:net.nₜ+net.nₚ])"
    end
    t, r, p, ψ
end

"""
Simulate a network and return the full time series.
- i: input filename, e.g. "net.bson"
- mut_id: index of a node to mutate, i.e. knock out.
- o: fname. output files with time along axis 2 and protein along axis 1
- t0: fname. starting values. Use this to let a mutated network start where the wildtype converged. 
    Columns named r_NODENAME1, r_NODENAME2, ..., p_NODENAME1, ... ψ_NODENAME1, ...
    Node names are assumed to be in order as the nodes in the input network.
    Make them with either another simulate call, a steady state call or write them manually.
"""
@main function timeseries(i::String, mut_id=nothing; o=nothing, duration=nothing, t0=nothing)
    # default output name
    if o === nothing
        o = mut_id === nothing ? "sim.tsv" : "sim_mut$mut_id.tsv"
    end

    net = loadnet(i)
    nₜₚ = net.nₜ+net.nₚ
    
	if t0 === nothing
        r0, p0, ψ0 = nothing, nothing, nothing
    else
        _, r0, p0, ψ0 = read_timeseries(t0, net)
        # final values
        r0, p0, ψ0 = r0[end, :], p0[end, :], ψ0[end, :]
    end
	
	u₀ = ODEs.get_u₀(net, r0, p0, ψ0)
    @assert !any(isnan.(u₀))
	solution = ODEs.@domainerror ODEs.timeseries(net, mut_id, u₀, duration)
	solution === nothing && return
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
        # dimensions are (node, RNA/protein/activation, time)
        # transpose to have time along first dim.
        r = solution[:, 1, :]'
        p = solution[:, 2, :]'
		ψ = solution[1:nₜₚ, 3, :]'
        t = reshape(solution.t, :, 1) # make it column vector here
        savedlm(o, hcat(t, r, p, ψ); colnames=[["t"]; "r_" .* net.names; "p_" .* net.names; "ψ_" .* net.names[1:nₜₚ]])
	end
end

"""
Plot time series.
- i: fname e.g. "sim.tsv". size = #times × (mRNA + protein + phos measurements)
- nₜ: number of TFs.
- nₚ: number of KPs.
- o: optional file to write plot to
"""
@main function plot(i::String, nₜ::Integer, nₚ::Integer; o=stdout)
    include(SRC * "simulation/PlotTimeSeries.jl") # save startup time by only including it when necessary
    t, r, p, ψ = read_timeseries(i)
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
- i: network fname, e.g. "net.bson"
- mut_id: index indicating which mutation to perform.
If "mutations" is not provided, it refers to index of the protein to mutate.
- r: out fname for r values.
- p: out fname for p values.
- psi: out fname for ψ values.
- mut_file: optional fname to provided where "mut_id" will then be the index of a column.
If "mut" is not provided, the first (and ideally only) column of the file will be used.
"""
@main function steadystate(i::String, mut_id=nothing; r=nothing, p=nothing, psi=nothing, mut_file=nothing)
    # default output names
    mutSuf = mut_id === nothing ? "" : "_mut$mut_id"
	if r === nothing r = "steady_r$mutSuf.mat" end
	if p === nothing p = "steady_p$mutSuf.mat" end
	if psi === nothing psi = "steady_psi$mutSuf.mat" end

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
- net: a network filename, e.g. "net.bson"
- o: stdout or file to write result to
"""
@main function logFC(net::String; o=stdout)
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

