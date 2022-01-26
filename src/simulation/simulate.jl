#!/usr/bin/env julia
@src "simulation/GeneRegulation"
@src "simulation/ODEs"
@use "utilities/ReadWrite"
using Fire
using LinearAlgebra
using Test  # @test_logs used in randomNetLogFCs
using .Threads: @threads
using DataFrames
using Plots: savefig

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
    else
        @assert all(chop.(colnames[r0idx], head=2, tail=0) .== chop.(colnames[p0idx], head=2, tail=0)) "colnames does not match between r_COLNAME and p_COLNAME"
        @assert all(chop.(colnames[p0idx], head=2, tail=0)[1:sum(ψ0idx)] .== chop.(colnames[ψ0idx], head=2, tail=0)) "colnames does not match between p_COLNAME and ψ_COLNAME"
    end
    node_names = chop.(colnames[r0idx], head=2, tail=0)
    t, r, p, ψ, node_names
end

function read_u₀(fname::String, net=nothing)
    _, r0, p0, ψ0, _ = read_timeseries(fname, net)
    # final values
    r0, p0, ψ0 = r0[end, :], p0[end, :], ψ0[end, :]
	u₀ = ODEs.get_u₀(net, r0, p0, ψ0)
    @assert !any(isnan.(u₀))
    u₀
end
function read_u₀(::Nothing, net=nothing)
	u₀ = ODEs.get_u₀(net)
    @assert !any(isnan.(u₀))
    u₀
end


"""
Simulate a network and return the full time series.
- i: input filename, e.g. "net.bson"
- mut_id: index of a node to mutate, i.e. knock out.
- o: fname. output files with time along axis 2 and protein along axis 1
- u0: fname. starting values. Use this to let a mutated network start where the wildtype converged. 
    Columns named r_NODENAME1, r_NODENAME2, ..., p_NODENAME1, ... ψ_NODENAME1, ...
    Node names are assumed to be in order as the nodes in the input network.
    Make them with either another simulate call, a steady state call or write them manually.
"""
@main function timeseries(i::String, mut_id=nothing; o=nothing, duration=nothing, u0=nothing)
    # default output name
    if o === nothing
        o = mut_id === nothing ? "sim.tsv" : "sim_mut$mut_id.tsv"
    end

    net = loadnet(i)
    nₜₚ = net.nₜ+net.nₚ
    u₀ = read_u₀(u0, net)

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

# create a color hex that is the same as "initial" except two characters are randomly changed.
function color_variant(initial::String)
    valid = ['0':'9'; 'a':'f']
    chars = collect(initial)
    for replace_ind in [3, 5, 7]
        chars[replace_ind] = rand(setdiff(valid, chars[replace_ind]))
    end
    join(chars)
end

"""
Plot time series.
- i: fname e.g. "sim.tsv". size = #times × (mRNA + protein + phos measurements)
- nₜ: number of TFs.
- nₚ: number of KPs.
- o: optional file to write plot to
"""
@main function plot_timeseries(i::String, nₜ::Integer, nₚ::Integer; o=stdout)
    @src "simulation/PlotTimeSeries"
    t, r, p, ψ, gene_names = read_timeseries(i)
    # dimensions should be time points along first dim, nodes along second in all files.

    dfs = []
    for (i, name) ∈ enumerate(gene_names)
        type = i <= nₜ ? "TF" : (i <= nₜ+nₚ ? "KP" : "O")
        push!(dfs, DataFrame("value"=>r[:, i], "time"=>t, "type"=>type, "name"=>name, "measure"=>"mRNA"))
        push!(dfs, DataFrame("value"=>p[:, i], "time"=>t, "type"=>type, "name"=>name, "measure"=>"prot"))
        i > nₜ+nₚ || push!(dfs, DataFrame("value"=>ψ[:, i], "time"=>t, "type"=>type, "name"=>name, "measure"=>"phos"))
    end
    for i in 1:length(dfs) dfs[i][!, "series"] .= i end

    df = vcat(dfs...)
    df[!, "label"] = df[!, "name"] .* " " .* df[!, "measure"]
    df[!, "style"] .= :solid; df[df[!, "measure"] .== "phos", "style"] .= :dash
    df[!, "width"] .= 1; df[df[!, "measure"] .== "mRNA", "width"] .= 3
    df[!, "subplot"] = df[!, "type"]
    df[!, "color"] .= "#000000"
    # choose colors
    type2color = Dict("KP"=>"#B663B1", "TF"=>"#6DB845", "O"=>"#545000")
    colors_TF = Dict(n=>type2color["TF"] for n ∈ gene_names[1:nₜ])
    colors_KP = Dict(n=>type2color["KP"] for n ∈ gene_names[nₜ+1:nₜ+nₚ])
    colors_O  = Dict(n=>type2color["O"]  for n ∈ gene_names[nₜ+nₚ+1:end])
    name2color = Dict(colors_TF..., colors_KP..., colors_O...)
    # add variation
    name2color = Dict(n=>color_variant(c) for (n, c) ∈ name2color)
    # add to table
    df[!, "color"] .= [name2color[n] for n ∈ df[!, "name"]]

	p = PlotTimeSeries.plot_timeseries(df)

	if o != stdout
        savefig(p, o)
    else
        display(p)
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
- o: out tab-separated file with columns r, p, ψ and node names as rownames. Default "steady.tsv" or "steady_mut#.tsv".
- mut_file: optional fname to provided where "mut_id" will then be the index of a column.
If "mut" is not provided, the first (and ideally only) column of the file will be used.
"""
@main function steadystate(i::String, mut_id=nothing; o=nothing, mut_file=nothing, u0=nothing)
    # default output names
	if o === nothing
	    o = mut_id === nothing ? "steady.tsv" : "steady_mut$mut_id.tsv"
	end

	net = loadnet(i)
	solution = steady_state(net, mut_id, mut_file)
	@info(solution.retcode)
	if solution.retcode in [:Success, :Terminated]
		r = solution[:, 1, end]
		p = solution[:, 2, end]
		ψ = solution[:, 3, end] # including entries for O that are always zero so the written table will be rectangular.
        savedlm(o, hcat(r, p, ψ); colnames=["r", "p", "ψ"], rownames=net.names)
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
- o: file to write result to
"""
@main function logFC(net::String, o::String)
    net = loadnet(net)
	measurements = ODEs.@domainerror(ODEs.logFC(net))
	if measurements !== nothing
		@info("logFC values simulated")
        savedlm(o, measurements; colnames=net.names[1:net.nₜ+net.nₚ], rownames=net.names)
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
    incl("src/simulation/weight");  # random(...) and correct()
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

