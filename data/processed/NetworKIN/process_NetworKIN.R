
library(data.table)

setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/NetworKIN")

biogridPTM = fread("networkin_biogrid.tsv.gz", sep="\t", quote="", header=T, select=c(1,5:8), col.names=c("Target", "KP", "NetworKIN", "NetPhorest", "STRING"))
PTM = fread("networkin.tsv.gz", sep="\t", quote="", header=T, select=c(1,4:7), col.names=c("Target", "KP", "NetworKIN", "NetPhorest", "STRING"))

biogrid_edges = biogridPTM[,list(NetworKIN=max(NetworKIN), NetPhorest=max(NetPhorest), STRING=max(STRING)), by=list(KP, Target)]
edges = PTM[,list(NetworKIN=max(NetworKIN), NetPhorest=max(NetPhorest), STRING=max(STRING)), by=list(KP, Target)]

mapping = match(paste(biogrid_edges$KP, biogrid_edges$Target), paste(edges$KP, edges$Target))
cor(biogrid_edges$NetworKIN, edges$NetworKIN[mapping])
cor(biogrid_edges$NetPhorest, edges$NetPhorest[mapping])
cor(biogrid_edges$STRING, edges$STRING[mapping])
all(biogrid_edges$NetworKIN <= edges$NetworKIN[mapping])
all(biogrid_edges$NetPhorest <= edges$NetPhorest[mapping])
all(biogrid_edges$STRING == edges$STRING[mapping])

# all scores in the more untrusted larger set "edges" are larger than in the smaller more trusted set "edges_biogrid" 
# which was made by only predicting kinase binding to known phosphorylation sites from biogrid data

write.table(biogrid_edges, "scores_biogrid.tsv", sep="\t", row.names=F, quote=F)
write.table(edges, "scores.tsv", sep="\t", row.names=F, quote=F)






