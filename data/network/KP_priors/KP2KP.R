#!/usr/bin/env Rscript
# packages
library(data.table)
library(fdrtool)

setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")

# made in parts, takes half an hour:
KP2KP = fread("KP2KP_medianest.tsv", sep="\t", header=T)

KP2KP[,median_weight:= -estimate * sign]
KP2KP[,sign:= sign(median_weight)]
KP2KP[,q:=fdrtool(p.value, statistic="pvalue", plot=FALSE)$qval]
stopifnot(all(unique(KP2KP$KP) == KPs))

colnames(KP2KP)[colnames(KP2KP)=="KP_reg"] = "substrate"
KP2KP = KP2KP[,c("KP", "substrate", "p.value", "q", "sign", "median_weight")] # reorder and select columns
write.table(KP2KP, "KP2KP.tsv", sep="\t", quote=F, row.names=F)


# analysis
sum(KP2KP$q < .2)
