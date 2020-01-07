
options(stringsAsFactors=FALSE)
setwd("~/cwd/data/nodes")

# read TF edge datasets
nodes = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2]
nodes = rbind(nodes, read.table("../processed/balaji_2008/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T))
nodes = rbind(nodes, read.table("../processed/yeastract/TF_edges.tsv", sep="\t", quote="", check.names=F, header=F, col.names=c("TF", "Target", "Mode"))[,1:2])
# harbison has 3 kinases that has been analysed as if they are TFs. We simply make sure to let TF be reduced by excluding anything in the KP set
nodes = rbind(nodes, read.table("../processed/harbison_2004/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T, col.names=c("TF", "Target", "Pval"))[,1:2])
nodes = rbind(nodes, read.table("../processed/horak_2002/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2])
nodes = rbind(nodes, read.table("../processed/workman_2006/TF_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2])
TF = unique(nodes$TF)

# read KP edge datasets
nodes = read.table("../processed/biogrid/P_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2]
nodes = rbind(nodes, read.table("../processed/fasolo_2011/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T))
nodes = rbind(nodes, read.table("../processed/fiedler_2009/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T, col.names=c("P", "Target", "Mode"))[,1:2])
nodes = rbind(nodes, read.table("../processed/parca_2019/P_edges.tsv", sep="\t", quote="", check.names=F, header=T))
KP_edges = unique(nodes$P)

# read from perturbation files
perturbation1 = read.table("../processed/chua_2006/TF_KO.tsv", quote="", check.names=F, header=T, row.names=1)
perturbation2 = read.table("../processed/chua_2006/TF_OE.tsv", quote="", check.names=F, header=T, row.names=1)
perturbation3 = read.table("../processed/fiedler_2009/avg_KO.tsv", quote="", check.names=F, header=T, row.names=1)
perturbation4 = read.table("../processed/goncalves_2017/KP_KO_TF.tsv", quote="", check.names=F, header=T, row.names=1)
perturbation5 = read.table("../processed/goncalves_2017/KP_KO_KP.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
perturbation6 = read.table("../processed/holstege_2010/PK_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
perturbation7 = read.table("../processed/holstege_2014/KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
perturbation8 = read.table("../processed/luscombe_2010/TF_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
perturbation9 = read.table("../processed/zelezniak_2018/PK_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
perturbation10 = read.table("../processed/luscombe_2010/PK_KO.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1)
KP_KO = unique(c(colnames(perturbation3), colnames(perturbation4), colnames(perturbation5), colnames(perturbation9), colnames(perturbation10), unlist(strsplit(colnames(perturbation6), "_"))))
KP_meas = rownames(perturbation5)
KP = unique(c(KP_edges, KP_KO, KP_meas))
KP_KO = unique(c(KP_KO, KP[KP %in% colnames(perturbation7)])) # add the KP edges that KOed in Holstege 
KP = unique(c(KP_edges, KP_KO, KP_meas))

TF = c(TF, colnames(perturbation1), colnames(perturbation2), colnames(perturbation4), colnames(perturbation8))
TF = c(TF, rownames(perturbation4)) # not sure about the rownames added
measured = c(rownames(perturbation1), rownames(perturbation2), rownames(perturbation3), rownames(perturbation4), 
      rownames(perturbation6), rownames(perturbation7), rownames(perturbation8))
all(colnames(perturbation9) %in% measured)
all(rownames(perturbation9) %in% measured)

KP = sort(unique(KP)); KP_KO = sort(unique(KP_KO))
TF = sort(unique(TF[!(TF%in%KP)]))
measured = sort(unique(measured))

write.table(TF, file="TF_putative.txt", quote=F, row.names=F, col.names=F)
# a TF is only assigned to the set of TFs if it has a regulatory role, otherwise it can be in O but only if it is measured
edges = read.table("../TF_edges/edges.tsv", sep="\t", header=T, quote="", stringsAsFactors=F)
edges = edges[edges$Target %in% measured,]
TF = TF[TF %in% edges$TF]
O = measured[!(measured %in% c(KP,TF))]

write.table(  KP,      "KP.txt", quote=F, row.names=F, col.names=F)
write.table(  KP_KO,"KP_KO.txt", quote=F, row.names=F, col.names=F)
write.table(     TF,   "TF.txt", quote=F, row.names=F, col.names=F)
write.table(c(KP,TF,O), "V.txt", quote=F, row.names=F, col.names=F)



