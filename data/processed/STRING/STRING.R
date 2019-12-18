
library(reshape2)

setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/STRING")

directed_melt = read.table("directed.tsv", sep="\t", header=T, quote="")
edges = unique(directed_melt[,c("Source","Target")])
directed = edges
modes = list()
for (lvl in levels(directed_melt$Mode)) {
    modes[[lvl]] = directed_melt[directed_melt$Mode == lvl,]
    modes[[lvl]] = modes[[lvl]][,!(colnames(modes[[lvl]]) %in% c("Mode", "Action"))]
    modes[[lvl]] = aggregate(Score ~ Source + Target, data=modes[[lvl]], mean)
    colnames(modes[[lvl]])[colnames(modes[[lvl]]) == "Score"] = lvl
    directed = merge(directed, modes[[lvl]], all=T)
}

write.table(directed, "scores.tsv", sep="\t", quote=F, row.names=F, na="")
