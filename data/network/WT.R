#!/usr/bin/env Rscript
library(Matrix)
library(ggplot2)
library(fdrtool)
library(data.table)

options(stringsAsFactors=FALSE)

# functions
flatten = function(x) as.vector(as.matrix(x))
qhalfnorm_ = function(pvals, sd) {
    vals = qhalfnorm(edges$Pval, theta=sd2theta(sd), lower.tail=F)
    vals[vals == Inf] = max(vals[vals != Inf])
    vals
}
sparsematrix = function(i, j, x) {
    i = match(i, Vs)
    j = match(j, TFs$TF)
    dims_ = list(length(Vs), nrow(TFs))
    dimnames_ = list(Vs, TFs$TF)
    as.matrix(sparseMatrix(i=i, j=j, x=x, dims=dims_, dimnames=dimnames_))
}

add_noise = function(adjacency, noise_sd = 1/sqrt(prod(dim(adjacency)))) {
    cat(noise_sd, "\n")
    adjacency_noise = adjacency
    lacking = adjacency_noise==0
    adjacency_noise[lacking] = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))[lacking]
    adjacency_noise
}

get_sign_mask = function(adjacency) {
    out = sign(adjacency)
    out[out == +1] = "+"
    out[out == -1] = "-"
    out[out ==  0] = "."
    out
}

# read
edges = read.table("TF_priors/TF_edges.tsv", sep="\t", header=T, quote="")
TFs = read.table("TF_mode.tsv", sep="\t", header=T, quote="")
Vs  = flatten(read.table("V_protein.txt", quote=""))
nTF = nrow(TFs)
nV = length(Vs)

## get weight sign from mode of regulation. Some edges already have a mode. Most do not, so they will get their mode from the TF.
edge_modes = TFs$Mode[match(edges$TF, TFs$TF)]
edges$Mode[edges$Mode == ""] = edge_modes[edges$Mode == ""]
# write repressor instead of both repressor and inhibitor
edges$Mode[edges$Mode == "inhibitor"] = "repressor"
# sum(edges$Mode == "activator")
# sum(edges$Mode == "repressor")
edges$sign = 1; edges$sign[edges$Mode == "repressor"] = -1

# create edge weights
edges$gauss001 = qhalfnorm_(edges$Pval, 0.01) * edges$sign
edges$gauss01  = qhalfnorm_(edges$Pval, 0.10) * edges$sign
edges$gauss025 = qhalfnorm_(edges$Pval, 0.25) * edges$sign

edges_FDR10 = edges[edges$qval<.1,]
edges_FDR20 = edges[edges$qval<.2,]

### write

write.table(edges_FDR20, "TF_edge_weights.tsv", sep="\t", quote=F, row.names=F)

# Get as matrix: weights, p-values, sign, and weights with noise
get_adjacency = function(edges, weight) {
    sparsematrix(edges$Target, edges$TF, edges[,weight])
}

get_p_adjacency = function(edges, weight) {
    adjacency_pval = sparsematrix(edges$Target, edges$TF, edges$Pval)
    stopifnot(!any(edges$Pval==0))  # should be true then we can safely do:
    adjacency_pval[adjacency_pval == 0] = NaN  # these entries are only zero because they were that way by default, not because a p-value=0 has been copied here.
    adjacency_pval
}

WT_mask = get_sign_mask(get_adjacency(edges, "gauss001"))
write.table(WT_mask, "WT_mask.ssv", sep=" ", quote=F)

fwrite(sparsematrix(edges$Target, edges$TF, edges$gauss01), "WT_gauss01.ssv", sep=" ")
fwrite(sparsematrix(edges_FDR20$Target, edges_FDR20$TF, edges_FDR20$gauss001), "WT_FDR20_gauss001.ssv", sep=" ")


