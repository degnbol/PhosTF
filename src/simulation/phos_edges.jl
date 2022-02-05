#!/usr/bin/env julia
using Statistics: mean
using SparseArrays

# This file supplies init_Wₚ₊Wₚ₋ and uses random_λ. It is used in Network.jl and is just moved here to split up the work a bit.


"""
Convert a Wₚ from int to float.
Choosing the float values are done in a manner that will make it likely 
that presence or absence of regulation should make a difference to down-stream gene expression levels.
- Wₚ: Matrix of 1, -1, and 0
"""
function init_Wₚ₊Wₚ₋(genes::Vector, Wₚ::Matrix{<:Integer}, λ₊::Vector, λ₋::Vector; hyper::Dict=Dict())
    if length(hyper) == 0
        hyper = Dict("cancel"=>true, "mean_k"=>true, "vec"=>false, "norm"=>0)
    end

    nₚ = size(Wₚ,2)
    nₜ = size(Wₚ,1) - nₚ
    Wₚ = Float64.(Wₚ)
    
    Wₚ₊, Wₚ₋ = Wₚ₊Wₚ₋(Wₚ)
    
    # on each node there will be a decay that the edges will act against collectively.
    # In some cases it will be 0/0 hence the max(1,...)
    cancel_λ₊ = λ₊
    cancel_λ₋ = λ₋
    if hyper["cancel"]
        cancel_λ₊ ./= max.(1., sum(Wₚ₋ .> 0; dims=2))
        cancel_λ₋ ./= max.(1., sum(Wₚ₊ .> 0; dims=2))
    end
    
    # Each edge is given the value that would make sense if we imagined it was the only regulator of its target, 
    # which would be higher than the passive decay of a similar magnitude, so we use the same random distribution here.
    # Positive edges onto a node has to work against negative decay and vice versa.
      
    # KP->TF. Take into account how well each target TF binds to their regulon.
    μk = hyper["mean_k"] ? mean_k(genes, nₜ) : ones(nₜ)
    Wₚ₊[1:nₜ, :]       .*= cancel_λ₋[1:nₜ]       .+ (hyper["vec"] ? random_λ(nₜ) : random_λ(nₜ, nₚ)) .* μk
    Wₚ₋[1:nₜ, :]       .*= cancel_λ₊[1:nₜ]       .+ (hyper["vec"] ? random_λ(nₜ) : random_λ(nₜ, nₚ)) .* μk
    # KP->KP
    Wₚ₊[nₜ+1:nₜ+nₚ, :] .*= cancel_λ₋[nₜ+1:nₜ+nₚ] .+ (hyper["vec"] ? random_λ(nₚ) : random_λ(nₚ, nₚ))
    Wₚ₋[nₜ+1:nₜ+nₚ, :] .*= cancel_λ₊[nₜ+1:nₜ+nₚ] .+ (hyper["vec"] ? random_λ(nₚ) : random_λ(nₚ, nₚ))

    if hyper["norm"] > 0
        # Reduce effects so they are of similar total magnitude onto each target.
        # Uses max(1,x) to avoid nan div.
        if hyper["norm"] > 1
            Wₚ₋ ./= max.(1., sum(Wₚ₋; dims=2))
            Wₚ₊ ./= max.(1., sum(Wₚ₊; dims=2))
        else
            Wₚ₋ ./= max.(1., sum(Wₚ₋ .> 0; dims=2))
            Wₚ₊ ./= max.(1., sum(Wₚ₊ .> 0; dims=2))
        end
    end
    
    Wₚ₊, Wₚ₋
end


"""
Each gene is regulated by a number of modules. Each of those modules contains a number of inputs which are TFs.
Each of those regulations have a k value indicating the fraction of binding required for halv activation (dissociation constant k).
Looking through all these modules we can find any TF that regulates a gene, and we call it an edge. 
We find the k for each of these edges.
"""
function edge_ks(genes::Vector)
    n = length(genes)
    source = [i for g in 1:n for m in genes[g].modules for i in m.inputs]
    target = [g for g in 1:n for m in genes[g].modules for i in m.inputs]
    value  = [k for g in 1:n for m in genes[g].modules for k in m.k]
    sparse(target, source, value, n, n, mean)
end
"""
avg. dissociation constant k for the outgoing edges of each TF.
"""
function mean_k(genes::Vector, nₜ::Integer)
    ks = edge_ks(genes)[:, 1:nₜ]
    # average the nonzero entries in each row
    nk = sum(ks .!= 0, dims=1) |> vec
    ks = sum(ks, dims=1) |> vec
    ks[nk .> 0] ./= nk[nk .> 0]
    ks
end
