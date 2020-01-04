
library(reshape2)

setwd("~/cwd/data/processed/STRING")

directed_melt = read.table("directed.tsv", sep="\t", header=T, quote="")
edges = unique(directed_melt[,c("Source","Target")])
directed = edges
types = list()
for (lvl in levels(directed_melt$Mode)) {
    types[[lvl]] = directed_melt[directed_melt$Mode == lvl,]
    types[[lvl]] = types[[lvl]][,!(colnames(types[[lvl]]) %in% c("Mode", "Action"))]
    types[[lvl]] = aggregate(Score ~ Source + Target, data=types[[lvl]], mean)
    colnames(types[[lvl]])[colnames(types[[lvl]]) == "Score"] = lvl
    directed = merge(directed, types[[lvl]], all=T)
}

write.table(directed, "scores.tsv", sep="\t", quote=F, row.names=F, na="")

mode_scores = directed[,c("Source", "Target", "activation", "inhibition")]
mode_scores = mode_scores[!(is.na(mode_scores$activation) & is.na(mode_scores$inhibition)),]
mode_scores$Mode = NA
mode_scores$Mode[mode_scores$activation > mode_scores$inhibition | is.na(mode_scores$inhibition)] = "activator"
mode_scores$Mode[mode_scores$activation < mode_scores$inhibition | is.na(mode_scores$activation)] = "inhibitor"
mode_scores = mode_scores[,c(1,2,5)]
mode_scores = na.omit(mode_scores)

write.table(mode_scores, "modes.tsv", sep="\t", quote=F, row.names=F)

