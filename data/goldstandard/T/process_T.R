
setwd("~/cwd/data/goldstandard/T")

balaji = read.table("balaji.tsv", sep="\t", col.names=c("TF", "Target"))
balaji_2008 = read.table("balaji_2008.tsv", sep="\t", col.names=c("TF", "Target"))
harbison = read.table("harbison.tsv", sep="\t", col.names=c("TF", "Target"))
horak = read.table("horak.tsv", sep="\t", col.names=c("TF", "Target"))
workman = read.table("workman.tsv", sep="\t", header=T)
yeastract = read.table("yeastract.tsv", sep="\t", col.names=c("TF", "Target", "Mode"), na.strings="-")
STRING = read.table("../../processed/STRING/modes.tsv", sep="\t", header=T); colnames(STRING)[1] = "TF"

edges = unique(rbind(balaji, balaji_2008, harbison, horak, workman, yeastract[,1:2]))
# add regulation mode
edges = merge(edges, yeastract, all.x=T)
edges = merge(edges, STRING, all.x=T)

write.table(edges, file="edges.tsv", sep="\t", row.names=F, quote=F, na="-")


