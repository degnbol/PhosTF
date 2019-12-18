
library(reshape2)

setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/harbison_2004")

YPD = as.matrix(read.table("YPD.tsv", header=T, row.names=1, sep="\t", check.names=F))
conds = as.matrix(read.table("conds.tsv", header=T, row.names=1, sep="\t", check.names=F))

YPD_melt = melt(YPD); colnames(YPD_melt) = c("gene", "TF", "Pval")
conds_melt = melt(conds); colnames(conds_melt) = c("gene", "TF", "Pval")

pvals = aggregate(Pval ~ ., data=rbind(YPD_melt, conds_melt), min)

write.table(pvals[,c(2,1,3)], file="TF_edges.tsv", quote=F, row.names=F, sep="\t")

