#!/usr/bin/env Rscript
# functions
flatten = function(x) as.vector(as.matrix(x))


setwd("/Users/christian/GoogleDrev/PKTFX/data/goldstandard/P")

PK = flatten(read.table("PK.txt", quote="", stringsAsFactors=F))
PP = flatten(read.table("PP.txt", quote="", stringsAsFactors=F))
P_edges = read.table("edges.tsv", quote="", sep="\t", header=T, stringsAsFactors=F)

P_edges$Relationship[(P_edges$P %in% PK) & (P_edges$Relationship == "-")] = "kinase"
P_edges$Relationship[(P_edges$P %in% PP) & (P_edges$Relationship == "-")] = "phosphatase"

# some numbers for sanity check
sum(P_edges$Relationship == "kinase")
sum(P_edges$Relationship == "phosphatase")
sum(P_edges$Relationship == "-")
nrow(P_edges)

write.table(P_edges, file="edges_resign.tsv", sep="\t", quote=F, row.names=F)


