#!/usr/bin/env Rscript

# settings
options(stringsAsFactors=F)

# packages
library(data.table)
library(fdrtool)

setwd("~/cwd/data/network/KP_priors")
# read
load("wilcoxon.RData")
TF_edges = fread("../TF_edge_weights.tsv", sep="\t", header=T)
KP_pert_J = fread("../../perturbation/J_inner.csv", sep=",", header=T)
KP_pert_J = melt(KP_pert_J[, c("ORF", KPs), with=FALSE], id.vars="ORF")
colnames(KP_pert_J)[2:3] = c("KP", "j")
stopifnot(all(KP_pert_J$ORF == KP_pert$ORF))
stopifnot(all(KP_pert_J$KP == KP_pert$KP))
KP_pert = KP_pert[(KP_pert_J$j == 0),]
stopifnot(!any(is.na(KP_pert$M)))  # NaNs should be removed with U=1-J mask

KP2TF = fread("KP2TF.tsv", sep="\t", header=T)
KP2TF_sig = KP2TF[KP2TF$q < .01, c("KP", "TF")]
TF2V_sig = TF_edges[TF_edges$qval < .01, c("TF", "Target")]
colnames(TF2V_sig) = c("TF", "gene")  # rename to gene to match KP_pert

# collect all pathways
pathways = data.table()
for (kp in KPs) {
    tfs = KP2TF_sig[KP==kp,TF]
    pathways_ = TF2V_sig[TF%in%KP2TF_sig[KP==kp,TF],]
    pathways = rbind(pathways, data.table(gene=pathways_$gene, KP=kp))
}

# ignoring TF we reduce to unique KP->?->gene pathways
pathways = unique(pathways)
pathways[,Regulon:=TRUE]

KP_regulons = data.table(expand.grid(KP=KPs, gene=genes))
KP_regulons = KP_regulons[as.character(KP) != as.character(gene),]
KP_regulons[pathways, on=c("KP", "gene"), Regulon:=Regulon]
KP_regulons$Regulon[is.na(KP_regulons$Regulon)] = FALSE
colnames(KP_regulons)[colnames(KP_regulons)=="KP"] = "KP_reg"


kps = commandArgs(trailingOnly=TRUE)


for (kp in kps) {
    cat(kp, "\n")
    
    DT = KP_regulons # copy
    DT[KP_pert[KP== kp,], on="gene", M:=M]
    DT = na.omit(DT)
    # invert Regulon bool since case should be tested for being greater than control, and case is 0 while control is 1 apparently
    # it's very clear which way to test, the values reveal nothing significant if chosen in the wrong direction.
    DT = DT[, c(wilcox.test(abs(M) ~ !Regulon, alternative="g", conf.int=TRUE)[c("p.value", "estimate")], 
                .(sign=sign(wilcox.test(M ~ !Regulon, conf.int=TRUE)$estimate))), by=KP_reg]
    
    DT$KP = kp
    write.table(DT, paste0("KP2KP_parts/", kp, ".tsv"), sep="\t", quote=F, row.names=F)
}








