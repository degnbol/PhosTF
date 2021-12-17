#!/usr/bin/env julia
isdefined(Main, :Model) || include("../inference/Model.jl") 
isdefined(Main, :ReadWrite) || include("../utilities/ReadWrite.jl")
isdefined(Main, :ArrayUtils) || include("../utilities/ArrayUtils.jl")
isdefined(Main, :GraphUtils) || include("../utilities/GraphUtils.jl")


using Fire

# defaults
default_Wₜ, default_Wₚ = "WT.mat", "WP.mat"
default_net = "net.bson"


"Load file(s) as a single 2D array regardless if they match in length along axis 1."
hcatpad_load(fnames::Vector) = ArrayUtils.hcatpad(loaddlm(fname, Float64) for fname in fnames)
hcatpad_load(fname::String) = ReadWrite.loaddlm(fname, Float64)

"""
Write a graph defined by weight matrices to xgmml format.
- X: node values. Each column of X is used for a separate copy of the graph.
"""
@main function xgmml(Wₜ, Wₚ::String; o=stdout, title=nothing, X=[])
	Wₜ, Wₚ = ReadWrite.loaddlm(Wₜ), ReadWrite.loaddlm(Wₚ)
    nₜ, nₚ = size(Wₜ, 2), size(Wₚ, 2)
	xgmml((Wₜ, Wₚ), o, nₜ, nₚ, title, X)
end
@main function xgmml(i=default_net; o=stdout, title=nothing, X=[])
	net = ReadWrite.load(i, Network)
	xgmml(net, o, net.nₜ, net.nₚ, title, X)
end

"""
- net: either (Wₜ, Wₚ) or Network.
"""
function xgmml(net, o, nₜ, nₚ, title=nothing, X=[])
    title !== nothing || (title = o == stdout ? "PhosTF" : splitext(basename(o))[1])
	if isempty(X) write(o, GraphUtils.xgmml(net; title=title))
	else
		X = hcatpad_load(X)
		K = size(X, 2)
		# highlight each of the proteins if there are as many experiments as TFs+KPs
		highlight = nₜ+nₚ == K ? (1:K) : nothing
		write(o, GraphUtils.xgmml(net, X, highlight; title=title))
	end
end

