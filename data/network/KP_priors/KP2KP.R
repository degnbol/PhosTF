
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

KP2TF = fread("KP2TF_p.tsv", sep="\t", header=T)
KP2TF_sig = KP2TF[KP2TF$q < .1, c("KP", "TF")]
TF2V_sig = TF_edges[TF_edges$qval < .1, c("TF", "Target")]
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
pathways[KP_pert, on=c("KP", "gene"), M:=M]
pathways = na.omit(pathways)

KP2KP = data.table()
for (kp_substrate in KPs) {
    pathways_ = pathways  # copy
    # i. means we use the M from i ( [i, j, ...] )
    pathways_[pathways_[KP==kp_substrate], on="gene", substrate_M:=i.M]
    # na are rows for genes that kp_substrate does not regulate
    pathways_ = pathways_[!is.na(substrate_M), cor.test(M,substrate_M)[c("p.value", "estimate")], by=KP]
    pathways_[,substrate:=kp_substrate]
    KP2KP = rbind(KP2KP, pathways_)
}
KP2KP = KP2KP[KP!=substrate,]
KP2KP[,sign:=sign(estimate)]
KP2KP[,q:=fdrtool(p.value, statistic="pvalue", plot=FALSE)$qval]

KP2KP = KP2KP[,c("KP", "substrate", "estimate", "sign", "p.value", "q")]
write.table(KP2KP, "KP2KP.tsv", sep="\t", row.names=F, quote=F)

# analysis
hist(KP2KP[q<.2 & abs(estimate) > .2, .N, by=KP]$N)

