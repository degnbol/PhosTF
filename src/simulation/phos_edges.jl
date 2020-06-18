using Statistics
using SparseArrays

"""
Convert a Wₚ from int to float.
Choosing the float values are done in a manner that will make it likely 
that presence or absence of regulation should make a difference to down-stream gene expression levels.
"""
function init_Wₚ₊Wₚ₋(genes::Vector, Wₚ::Matrix, λ_phos::Vector)
    nₚ = size(Wₚ,1); nₜ = size(Wₚ,2) - nₚ
    Wₚ = 1.0Wₚ # convert type
    
    # here each edge is given the value that would make sense if we imagined it was the only regulator of its target
    Wₚ[1:nₚ,:] .*= 2λ_phos[1:nₚ,:]
    # mean_k = avg. dissociation constant k for the outgoing edges of each TF
    Wₚ[nₚ+1:nₚ+nₜ,:] .*= λ_phos[nₚ+1:nₚ+nₜ,:] .* mean_k(genes, nₜ, nₚ)
    
    Wₚ₊, Wₚ₋ = Wₚ₊Wₚ₋(Wₚ)

    # split the total effect evenly among edges that share a target. uses max(1,x) to avoid nan div
    Wₚ₋ ./= max.(1., sum(Wₚ₋ .> 0, dims=2))
    Wₚ₊ ./= max.(1., sum(Wₚ₊ .> 0, dims=2))
    
    Wₚ₊, Wₚ₋
end

function mean_edge_k(genes::Vector)
    n = length(genes)
    source = [i for g in 1:n for m in genes[g].modules for i in m.inputs]
    target = [g for g in 1:n for m in genes[g].modules for i in m.inputs]
    value  = [k for g in 1:n for m in genes[g].modules for k in m.k]
    sparse(target, source, value, n, n, mean)
end
function mean_k(genes::Vector, nₜ::Integer, nₚ::Integer)
    ks = mean_edge_k(genes)[:,nₚ+1:nₚ+nₜ]
    # average the nonzero entries in each row
    nk = sum(ks .!= 0, dims=1) |> vec
    ks = sum(ks, dims=1) |> vec
    ks[nk .> 0] ./= nk[nk .> 0]
    ks
end
