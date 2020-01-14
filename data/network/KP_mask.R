#!/usr/bin/env Rscript
# make a KP mask by limiting which KPs and TFs can be a target from the documented phosphorylation sites from BioGRID

# functions
flatten = function(x) as.vector(as.matrix(x))
read.flat = function(x) flatten(read.table(x))

setwd("~/cwd/data/network")

KPs = read.flat("KP.txt")
TFs = read.flat("TF.txt")
KPTFs = c(KPs, TFs)
targets = read.flat("../processed/biogrid/targets.txt")

mask = matrix(0, nrow=length(KPTFs), ncol=length(KPs))
mask[KPTFs %in% targets,] = 1
dimnames(mask) = list(KPTFs, KPs)

write.table(mask, "KP_mask.csv", sep=",", quote=F)
write.table(mask, "KP_mask.mat", sep=" ", quote=F, col.names=F, row.names=F)
