
options(stringsAsFactors=FALSE)
setwd("/Users/christian/GoogleDrev/PKTFX/data/nodes")

# read TF edge datasets
nodes = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T)
nodes = rbind(nodes, read.table("../processed/balaji_2008/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T))
nodes = rbind(nodes, read.table("../processed/yeastract/TF_edges.tsv", sep="\t", quote="", check.names=F, header=F, col.names=c("TF", "Target", "Mode"))[,1:2])
# harbison has 3 kinases that has been analysed as if they are TFs. We simply make sure to let TF be reduced by excluding anything in the KP set
nodes = rbind(nodes, read.table("../processed/harbison_2004/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T, col.names=c("TF", "Target", "Pval"))[,1:2])
nodes = rbind(nodes, read.table("../processed/horak_2002/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2])
nodes = rbind(nodes, read.table("../processed/workman_2006/TF_edges_nodup.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2])

TF = unique(nodes$TF)
V  = unique(c(TF, nodes$Target))

# read KP edge datasets
nodes = read.table("../processed/biogrid/P_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2]
nodes = rbind(nodes, read.table("../processed/fasolo_2011/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T))
nodes = rbind(nodes, read.table("../processed/fiedler_2009/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T, col.names=c("P", "Target", "Mode"))[,1:2])
nodes = rbind(nodes, read.table("../processed/parca_2019/P_edges.tsv", sep="\t", quote="", check.names=F, header=T))

KP = unique(nodes$P)
V = unique(c(V, nodes$Target))

# read from perturbation files
perturbation = read.table("../processed/chua_2006/TF_KO.tsv", quote="", check.names=F, header=T, row.names=1)
TF = c(TF, colnames(perturbation))
V = c(V, rownames(perturbation))
perturbation = read.table("../processed/chua_2006/TF_OE.tsv", quote="", check.names=F, header=T, row.names=1)
TF = c(TF, colnames(perturbation))
V = c(V, rownames(perturbation))
perturbation = read.table("../processed/fiedler_2009/avg_KO.tsv", quote="", check.names=F, header=T, row.names=1)[,2:3]
all(colnames(perturbation) %in% KP)
V = c(V, rownames(perturbation))
perturbation = read.table("../processed/goncalves_2017/KP_KO_TF.tsv", quote="", check.names=F, header=T, row.names=1)
KP = c(KP, colnames(perturbation))
# not sure about this one
TF = c(TF, rownames(perturbation))

perturbation = read.table("../processed/goncalves_2017/KP_KO_KP.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
KP = c(KP, colnames(perturbation), rownames(perturbation))
perturbation = read.table("../processed/holstege_2010/PK_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
KP = c(KP, unlist(strsplit(colnames(perturbation), "_")))
V = c(V, rownames(perturbation))
perturbation = read.table("../processed/luscombe_2010/TF_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
TF = c(TF, colnames(perturbation))
V = c(V, rownames(perturbation))
perturbation = read.table("../processed/zelezniak_2018/PK_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
KP = c(KP, colnames(perturbation))
all(colnames(perturbation) %in% V)
all(rownames(perturbation) %in% V)

KP = sort(unique(KP))
TF = sort(unique(TF[!(TF%in%KP)]))
V = sort(unique(c(KP,TF,V)))

write.table(KP, file="KP.txt", quote=F, row.names=F, col.names=F)
write.table(TF, file="TF.txt", quote=F, row.names=F, col.names=F)
write.table(V,  file="V.txt",  quote=F, row.names=F, col.names=F)





