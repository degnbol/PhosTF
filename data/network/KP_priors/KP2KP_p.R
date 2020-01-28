
# settings
options(stringsAsFactors=F)

# packages
library(data.table)

setwd("~/cwd/data/network/KP_priors")
# read
load("wilcoxon.RData")
KP_pert_J = fread("../../perturbation/J_inner.csv", sep=",", header=T)
KP_pert_J = melt(KP_pert_J[, c("ORF", KPs), with=FALSE], id.vars="ORF")
colnames(KP_pert_J)[2:3] = c("KP", "j")
stopifnot(all(KP_pert_J$ORF == KP_pert$ORF))
stopifnot(all(KP_pert_J$KP == KP_pert$KP))
KP_pert = KP_pert[(KP_pert_J$j == 0),]


KP2TF = read.table("KP2TF_p.tsv", sep="\t", quote="", header=T)
KP2TF_sig = KP2TF[KP2TF$p < 1e-2,]
TF2V_sig = TF_edges[TF_edges$Pval < 1e-2,]

KP2KP = expand.grid(Target=KPs, KP=KPs)[,2:1]
KP2KP$p=NA; KP2KP$sign=NA
KP2KP = KP2KP[KP2KP$KP != KP2KP$Target,]  # no self-loops
for (i in 1:nrow(KP2KP)) {
    source = KP2KP$KP[i]
    target = KP2KP$Target[i]
    regulon_TF = KP2TF_sig$TF[KP2TF_sig$KP == target]
    regulon    = unique(TF2V_sig$Target[TF2V_sig$TF %in% regulon_TF])
    if (length(regulon) > 0) {
        source_pert = KP_pert$M[KP_pert$KP == source & !(KP_pert$gene %in% regulon)]
        target_pert = KP_pert$M[KP_pert$KP == target & !(KP_pert$gene %in% regulon)]
        slope = summary(lm(source_pert ~ target_pert -1))$coefficients
        KP2KP$p[i] = slope[length(slope)]
        KP2KP$sign[i] = sign(slope[1])
    }
}


KP2KP$p_adj = p.adjust(KP2KP$p)

sum(KP2KP$p < 5e-2, na.rm=T)
sum(KP2KP$p < 1e-3, na.rm=T)
sum(KP2KP$p < 1e-6, na.rm=T)
sum(KP2KP$p < 1e-12, na.rm=T)
sum(KP2KP$p < 1e-15, na.rm=T)
sum(KP2KP$p < 1e-18, na.rm=T)
sum(KP2KP$p < 1e-21, na.rm=T)
hist(qbeta(na.omit(KP2KP$p), 1, 20, lower.tail=F), breaks=50)
hist(qbeta(na.omit(KP2KP$p[KP2KP$p < 0.05]), 1, 20, lower.tail=F), breaks=50)
sum(KP2KP$p_adj < 5e-2, na.rm=T)
sum(KP2KP$p_adj < 1e-6, na.rm=T)
hist(qbeta(na.omit(KP2KP$p_adj), 1, 20, lower.tail=F), breaks=50)
hist(qbeta(na.omit(KP2KP$p_adj[KP2KP$p_adj < 0.05]), 1, 20, lower.tail=F), breaks=50)



