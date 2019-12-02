
# functions
flatten = function(x) as.vector(as.matrix(x))


setwd("/Users/christian/GoogleDrev/PKTFX/data/goldstandard/P")

PK = flatten(read.table("PK.txt", quote="", stringsAsFactors=F))
PP = flatten(read.table("PP.txt", quote="", stringsAsFactors=F))
P_edges = read.table("edges.tsv", quote="", sep="\t", header=T, stringsAsFactors=F)
T_edges = read.table("../T/edges.tsv", quote="", sep="\t", header=T, stringsAsFactors=F)

# filter so nodes are only in P if they are not silent
n_edges = -1
# refilter, until filtering doesn't remove anything
P_edges_reduced = P_edges
while (n_edges != nrow(P_edges_reduced)) {
    n_edges = nrow(P_edges_reduced)
    # a P edge is valid if the target is a valid node in P or valid node in T
    regulate_P = P_edges_reduced$Target %in% P_edges_reduced$P
    regulate_T = P_edges_reduced$Target %in% T_edges$TF
    P_edges_reduced = P_edges_reduced[regulate_P | regulate_T,]
    cat("reduced to", nrow(P_edges_reduced), "\n")
}

write.table(P_edges_reduced, file="detectable/edges.tsv", sep="\t", quote=F, row.names=F)




