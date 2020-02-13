#!/usr/bin/env Rscript

# settings
options(stringsAsFactors=T)  # in order for merge to order according to KP TF O
# packages
library(data.table)
# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
perturbation = as.matrix(read.table("../../perturbation/logFC_inner.csv", sep=",", header=T, quote=""))
# all p-values are significant in TF_edges.tsv
TF_edges = read.table("../TF_edges.tsv", sep="\t", header=T, quote="")[,1:3] # discard incomplete edge mode column
TF_edges = TF_edges[TF_edges$TF %in% colnames(perturbation),]
KPs = flatten(read.table("../KP.txt"))
KPs = KPs[KPs %in% colnames(perturbation)]  # we can only do wilcoxon test for KPs that are pertubed
TFs = read.table("../TF_mode.tsv", sep="\t", header=T, quote="", stringsAsFactors=FALSE)
TFs = TFs[TFs$TF%in%TF_edges$TF,]
genes = rownames(perturbation)
nTF = nrow(TFs)
nKP = length(KPs)
N = length(genes)

KP_pert = melt(data.table(perturbation[,KPs], keep.rownames=T), id.vars="rn")
colnames(KP_pert) = c("gene", "KP", "M")

TF_regulons = data.table(expand.grid(TF=TFs$TF, Target=genes))
TF_regulons = merge(TF_regulons, TF_edges, all.x=T)
TF_regulons$Regulon = !is.na(TF_regulons$Pval)
TF_regulons[,Pval:=NULL]

# are they ordered KP, TF, O?
stopifnot(all(KP_pert$gene == rep(genes, nKP)))
stopifnot(all(TF_regulons$Target == rep(genes, nTF)))

write.table(KPs, "KP.txt", col.names=F, row.names=F, quote=F)
save.image("wilcoxon.RData")
