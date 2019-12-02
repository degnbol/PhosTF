setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/fiedler_2009")

edges = read.table("P_edges.tsv", header=T, sep="\t", quote="", check.names=F)
EMAP = read.table("EMAP_edges.tsv", header=T, sep="\t", quote="", check.names=F)

# collect both source, target and relationship since all have to match for the edge to be the same
edge_match = paste(edges$ORF, edges$`Target ORF`, edges$Relationship)
EMAP_match = paste(EMAP$ORF, EMAP$`Target ORF`, EMAP$Relationship)

edges$EMAP = EMAP$`EMAP Score`[match(edge_match, EMAP_match)]

write.table(edges, file="P_edges_EMAP.tsv", row.names=F, quote=F, sep="\t", na="NaN")
