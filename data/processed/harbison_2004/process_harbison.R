
library(reshape2)

setwd("~/cwd/data/processed/harbison_2004")

YPD = as.matrix(read.table("YPD.tsv", header=T, row.names=1, sep="\t", check.names=F))
YPD_lee_2002 = as.matrix(read.table("../lee_2002/YPD.tsv", header=T, row.names=1, sep="\t", check.names=F))
conds = as.matrix(read.table("conds.tsv", header=T, row.names=1, sep="\t", check.names=F))

YPD_melt = melt(YPD); colnames(YPD_melt) = c("gene", "TF", "Pval")
YPD_lee_2002_melt = melt(YPD_lee_2002); colnames(YPD_lee_2002_melt) = c("gene", "TF", "Pval")
conds_melt = melt(conds); colnames(conds_melt) = c("gene", "TF", "Pval")

pvals = aggregate(Pval ~ ., data=rbind(YPD_melt, YPD_lee_2002_melt, conds_melt), min)

write.table(pvals[,c(2,1,3)], file="TF_edges.tsv", quote=F, row.names=F, sep="\t")

