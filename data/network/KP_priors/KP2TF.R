#!/usr/bin/env Rscript
# packages
library(data.table)
library(fdrtool)

setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")

# made in parts, takes half an hour:
KP2TF = fread("KP2TF_medianest.tsv", sep="\t", header=T)

# the sign of medians is for the PK KO effect 
# so the sign in the file is KO_sign = - TF_sign * KP_sign => KP_sign = - KO_sign / TF_sign
# which for TF_sign \in {-1,1} is the same as KP_sign = - KO_sign * TF_sign
TFs$sign = (TFs$Mode == "activator") - (TFs$Mode == "repressor")

KP2TF[TFs[,c("TF","sign")], on="TF", TF_sign:=sign]
KP2TF[,median_weight:= - TF_sign * sign * estimate]
KP2TF[,sign:= sign(median_weight)]
KP2TF[,q:=fdrtool(p.value, statistic="pvalue", plot=FALSE)$qval]
stopifnot(all(unique(KP2TF$KP) == KPs))
KP2TF = KP2TF[,c("KP", "TF", "p.value", "q", "sign", "median_weight")] # reorder and select columns
write.table(KP2TF, "KP2TF.tsv", sep="\t", quote=F, row.names=F)


# analysis
sum(KP2TF$q < .2)
