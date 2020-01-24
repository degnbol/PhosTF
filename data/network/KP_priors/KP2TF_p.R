
# packages
library(data.table)

setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")

KP2TF = read.table("KP2TF_genesign.tsv", sep="\t", header=T)

# the sign is currently for the PK KO effect 
# so the sign in the file is KO_sign = - TF_sign * KP_sign => KP_sign = - KO_sign / TF_sign
# which for TF_sign \in {-1,1} is the same as KP_sign = - KO_sign * TF_sign
TFs$Mode[TFs$Mode == "activator"] = +1
TFs$Mode[TFs$Mode == "repressor"] = -1

KP2TF$sign = - TFs$Mode[match(KP2TF$TF, TFs$TF)] * KP2TF$sign

write.table(KP_edges, "KP2TF_p.tsv", sep="\t", quote=F, row.names=F)
