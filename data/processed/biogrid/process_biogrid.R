setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/biogrid")
# read
targets = read.table("BIOGRID-PTM-3.5.178.ptmtab.tsv", sep="\t", header=T, quote="")
sources = read.table("BIOGRID-PTM-RELATIONSHIPS-3.5.178.ptmrel.tsv", sep="\t", header=T, quote="")
mapping = match(sources$PTM.ID, targets$PTM.ID)
# combine
names(targets) = paste("Target", names(targets), sep=".")
edges = cbind(sources, targets[mapping,])
# remove NA
edges_nona = na.omit(edges)
cat((1- nrow(edges_nona) / nrow(edges)) * 100, "percent of edges were removed\n")
# sanity check
all(edges_nona$PTM.ID==edges_nona$Target.PTM.ID)
edges_nona = edges_nona[,names(edges_nona) != "Target.PTM.ID"]
names(edges_nona) = gsub("\\.", " ", names(edges_nona))
write.table(edges_nona, file="regulatory_edges.tsv", sep="\t", row.names=F, quote=F)
edges_core = edges_nona[edges_nona$`Target Post Translational Modification` == "Phosphorylation",]
edges_core = data.frame(P=edges_core$`Systematic Name`, Target=edges_core$`Target Systematic Name`, Relationship=edges_core$Relationship)
write.table(edges_core, file="P_edges.tsv", sep="\t", row.names=F, quote=F)
