
# packages
library(data.table)

setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")

# the sign is currently for the PK KO effect 
TFs$Mode[TFs$Mode == "activator"] = +1
TFs$Mode[TFs$Mode == "repressor"] = -1


KP_edges$sign[i] *= TFs$Mode[TFs$TF == TF]


write.table(KP_edges, "KP2TF_p.tsv", sep="\t", quote=F, row.names=F)
