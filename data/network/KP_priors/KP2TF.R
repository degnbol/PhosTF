#!/usr/bin/env Rscript
# packages
library(data.table)

setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")

# made in parts, takes half an hour:
KP2TF = read.table("KP2TF_medianest.tsv", sep="\t", header=T)

# the sign of medians is for the PK KO effect 
# so the sign in the file is KO_sign = - TF_sign * KP_sign => KP_sign = - KO_sign / TF_sign
# which for TF_sign \in {-1,1} is the same as KP_sign = - KO_sign * TF_sign
TFs$Mode[TFs$Mode == "activator"] = +1
TFs$Mode[TFs$Mode == "repressor"] = -1
TFs$Mode = as.numeric(TFs$Mode)

KP2TF$sign = - TFs$Mode[match(KP2TF$TF, TFs$TF)] * sign(KP2TF$estimate)
stopifnot(all(unique(KP2TF$KP) == KPs))
KP2TF$weight = - TFs$Mode[match(KP2TF$TF, TFs$TF)] * KP2TF$estimate
KP2TF = KP2TF[,c("KP", "TF", "p", "sign", "weight")] # reorder columns
write.table(KP2TF, "KP2TF.tsv", sep="\t", quote=F, row.names=F)
