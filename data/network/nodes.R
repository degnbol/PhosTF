#!/usr/bin/env Rscript
# define sets of KP, TF and O as well as filter TF edges by TF->V

options(stringsAsFactors=FALSE)
setwd("~/cwd/data/network")

# read TF edge datasets
TF_edges = read.table("TF_priors/TF_edges_putative.tsv", sep="\t", header=T, quote="", stringsAsFactors=F)
TF_edges_FDR = read.table("TF_priors/TF_edges_putative_FDR.tsv", sep="\t", header=T, quote="", stringsAsFactors=F)

# read KP edge datasets
KP_edges = read.table("../processed/biogrid/P_edges.tsv", sep="\t", quote="", check.names=F, header=T)[,1:2]
KP_edges = rbind(KP_edges, read.table("../processed/fasolo_2011/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T))
KP_edges = rbind(KP_edges, read.table("../processed/fiedler_2009/P_edges_ORF.tsv", sep="\t", quote="", check.names=F, header=T, col.names=c("P", "Target", "Mode"))[,1:2])
KP_edges = rbind(KP_edges, read.table("../processed/parca_2019/P_edges.tsv", sep="\t", quote="", check.names=F, header=T))
KP_with_edge = unique(KP_edges$P)

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
KP = unique(c(KP_with_edge, KP_KO, KP_meas))
KP_KO = unique(c(KP_KO, KP[KP %in% colnames(perturbation7)])) # add the KP edges that KOed in Holstege 
KP = unique(c(KP_with_edge, KP_KO, KP_meas))

measured = c(rownames(perturbation1), rownames(perturbation2), rownames(perturbation3), rownames(perturbation4), 
      rownames(perturbation6), rownames(perturbation7), rownames(perturbation8))
all(colnames(perturbation9) %in% measured)
all(rownames(perturbation9) %in% measured)

# a TF is only assigned to the set of TFs if it has a regulatory role, otherwise it can be in O but only if it is measured
TF = TF_edges$TF
# manually remove a few that are not genes but are promotor elements or protein complexes
TF = TF[!(TF %in% c("CSRE", "STRE", "ECB", "MCB", "SCB", "MSE", "SBF", "MATA2"))]

KP = sort(unique(KP))
KP_KO = sort(unique(KP_KO))
TF = sort(unique(TF[!(TF%in%KP)]))
measured = sort(unique(measured))

write.table(TF, file="TF_putative.txt", quote=F, row.names=F, col.names=F)
# filter TF edges. It is necessary to do (TF_edges$TF %in% TF) since KSS1 has been removed from TF by removal of KPs
TF_edges = TF_edges[(TF_edges$TF %in% TF) & (TF_edges$Target %in% c(KP,TF,measured)),]
TF_edges_FDR = TF_edges_FDR[(TF_edges_FDR$TF %in% TF) & (TF_edges_FDR$Target %in% c(KP,TF,measured)),]
TF_edges = TF_edges[TF_edges$TF != TF_edges$Target,] # no self-loops
TF_edges_FDR = TF_edges_FDR[TF_edges_FDR$TF != TF_edges_FDR$Target,] # no self-loops
stopifnot(all(TF %in% TF_edges$TF))  # if not true we should reduce set of TFs
TF[!TF%in%TF_edges_FDR$TF] # not all TFs are actually in use in the smaller FDR set of edges. Ignore for now but effictively some TFs are in V here.
hist(TF_edges$Pval)
hist(TF_edges_FDR$Pval)
write.table(TF_edges, file="TF_priors/TF_edges.tsv", sep="\t", row.names=F, quote=F, na="")
write.table(TF_edges_FDR, file="TF_priors/TF_edges_FDR.tsv", sep="\t", row.names=F, quote=F, na="")

# has to have TF edges onto, otherwise no way to regulate them
O = measured[!(measured %in% c(KP,TF)) & (measured %in% TF_edges$Target)]

write.table(  KP,      "KP.txt", quote=F, row.names=F, col.names=F)
write.table(  KP_KO,"KP_KO.txt", quote=F, row.names=F, col.names=F)
write.table(     TF,   "TF.txt", quote=F, row.names=F, col.names=F)
write.table(c(KP,TF,O), "V.txt", quote=F, row.names=F, col.names=F)

# add mode
KP_amigo = read.table("../processed/amigo2/KP.tsv", sep="\t", header=T)
TF_amigo = read.table("../processed/amigo2/TF.tsv", sep="\t", header=T)
KP = data.frame(KP=KP, Mode=KP_amigo$Mode[match(KP, KP_amigo$KP)])
TF = data.frame(TF=TF, Mode=TF_amigo$Mode[match(TF, TF_amigo$TF)])
KP_noKO = KP[!(KP$KP %in% KP_KO),]

write.table(KP, "KP.tsv", sep="\t", quote=F, row.names=F, na="")
write.table(TF, "TF.tsv", sep="\t", quote=F, row.names=F, na="")



