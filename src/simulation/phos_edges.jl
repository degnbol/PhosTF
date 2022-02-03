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
function init_Wₚ₊Wₚ₋(genes::Vector, Wₚ::Matrix{<:Integer}, λ₊::Vector, λ₋::Vector)
    nₚ = size(Wₚ,2)
    nₜ = size(Wₚ,1) - nₚ
    Wₚ = Float64.(Wₚ)
    
    Wₚ₊, Wₚ₋ = Wₚ₊Wₚ₋(Wₚ)

    # Each edge is given the value that would make sense if we imagined it was the only regulator of its target, 
    # which would be higher than the passive decay of a similar magnitude, so we use the same random distribution here.
    # Positive edges onto a node has to work against negative decay and vice versa.
      
    # KP->TF. Take into accont how well each target TF binds to their regulon.
    Wₚ₊[1:nₜ, :]       .*= λ₋[1:nₜ]       .+ random_λ(nₜ) .* mean_k(genes, nₜ)
    Wₚ₋[1:nₜ, :]       .*= λ₊[1:nₜ]       .+ random_λ(nₜ) .* mean_k(genes, nₜ)
    # KP->KP
    Wₚ₊[nₜ+1:nₜ+nₚ, :] .*= λ₋[nₜ+1:nₜ+nₚ] .+ random_λ(nₚ) 
    Wₚ₋[nₜ+1:nₜ+nₚ, :] .*= λ₊[nₜ+1:nₜ+nₚ] .+ random_λ(nₚ) 
    
    # Reduce effects so they are of similar total magnitude onto each target.
    # max(1,x) to avoid nan div.
    Wₚ₋ ./= max.(1., sum(Wₚ₋ .> 0, dims=2))
    Wₚ₊ ./= max.(1., sum(Wₚ₊ .> 0, dims=2))
    
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
