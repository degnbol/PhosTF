#!/usr/bin/env julia
isdefined(Main, :Model) || include("Model.jl") 
isdefined(Main, :ReadWrite) || include("src/utilities/ReadWrite.jl")
isdefined(Main, :ArrayUtils) || include("src/utilities/ArrayUtils.jl")
isdefined(Main, :Cytoscape) || include("src/Cytoscape.jl")


using Fire

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"


loadnet(i) = ReadWrite.load(i, Network)


"Load file(s) as a single 2D array regardless if they match in length along axis 1."
hcatpad_load(fnames::Vector) = ArrayUtils.hcatpad(loaddlm(fname, Float64) for fname in fnames)
hcatpad_load(fname::String) = ReadWrite.loaddlm(fname, Float64)

"""
Write a graph defined by weight matrices to xgmml format.
- X: node values. Each column of X is used for a separate copy of the graph.
"""
@main function xgmml(Wₜ, Wₚ::String; o=stdout, title=nothing, X=[])
	Wₜ, Wₚ = ReadWrite.loaddlm(Wₜ), ReadWrite.loaddlm(Wₚ)
	_,nₜ,nₚ = Model.nₒnₜnₚ(Wₜ,Wₚ)
	xgmml([Wₜ, Wₚ], o, nₜ, nₚ, title, X)
end
@main function xgmml(i=default_net; o=stdout, title=nothing, X=[])
	net = loadnet(i)
	xgmml([net], o, net.nₜ, net.nₚ, title, X)
end
@main function xgmml(W, nₜ::Integer, nₚ::Integer; o=stdout, title=nothing, X=[])
	W = ReadWrite.loaddlm(W)
	Wₜ, Wₚ = W[:,nₚ+1:nₚ+nₜ], W[:,1:nₚ] # we don't use the Model.WₜWₚ function since we want to allow P→X edges in case W==T
	xgmml([Wₜ, Wₚ], o, nₜ, nₚ, title, X)
end

function xgmml(i, o, nₜ, nₚ, title=nothing, X=[])
    title !== nothing || (title = o == stdout ? "pktfx" : splitext(basename(o))[1])
	if isempty(X) write(o, Cytoscape.xgmml(i...; title=title))
	else
		X = hcatpad_load(X)
		K = size(X,2)
		# highlight each of the proteins if there are as many experiments as PKs+TFs
		highlight = nₜ+nₚ == K ? (1:K) : nothing
		write(o, Cytoscape.xgmml(i..., X, highlight; title=title))
	end
end

