
library(reshape2)

setwd("~/cwd/data/processed/harbison_2004")

YPD = as.matrix(read.table("YPD.tsv", header=T, row.names=1, sep="\t", check.names=F))
YPD_lee = as.matrix(read.table("../lee_2002/YPD.tsv", header=T, row.names=1, sep="\t", check.names=F))
conds = as.matrix(read.table("conds.tsv", header=T, row.names=1, sep="\t", check.names=F))

# there are perfect overlap among the values in YPD from harbison 2004 and lee 2002 for the same experiments so we take the unique
YPD_melt = unique(na.omit(rbind(melt(YPD), melt(YPD_lee))))
conds_melt = na.omit(melt(conds))
colnames(YPD_melt) = c("Target", "TF", "Pval")
colnames(conds_melt) = c("Target", "TF", "Pval")

# number of measurements
nrow(YPD_melt)+nrow(conds_melt)

write.table(YPD_melt[,c(2,1,3)], "TF_edges_YPD.tsv", quote=F, row.names=F, sep="\t")
write.table(conds_melt[,c(2,1,3)], "TF_edges_conds.tsv", quote=F, row.names=F, sep="\t")
