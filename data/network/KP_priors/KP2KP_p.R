
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
stopifnot(!any(is.na(KP_pert$M)))  # NaNs should be removed with U=1-J mask

KP2TF = read.table("KP2TF_p.tsv", sep="\t", quote="", header=T)
KP2TF_sig = KP2TF[p.adjust(KP2TF$p, "BH") < 1e-3,]
TF2V_sig = TF_edges[p.adjust(TF_edges$Pval, "BH") < 1e-3,]

KP2KP = expand.grid(Target=KPs, KP=KPs)[,2:1]
KP2KP$p=NA; KP2KP$sign=NA
KP2KP = KP2KP[KP2KP$KP != KP2KP$Target,]  # no self-loops
for (i in 1:nrow(KP2KP)) {
    source = KP2KP$KP[i]
    target = KP2KP$Target[i]
    regulon_TF = KP2TF_sig$TF[KP2TF_sig$KP == target]
    regulon    = unique(TF2V_sig$Target[TF2V_sig$TF %in% regulon_TF])
    if (length(regulon) > 0) {
        source_pert = KP_pert[KP == source,]
        target_pert = KP_pert[KP == target,]
        # NaNs cause some KP->gene KO logFC to be missing
        nonaregulon = intersect(source_pert$gene, target_pert$gene)
        nonaregulon = nonaregulon[nonaregulon %in% regulon]
        source_pert = source_pert[gene%in%nonaregulon,M]
        target_pert = target_pert[gene%in%nonaregulon,M]
        # using correlation        
        if (length(nonaregulon) > 2) {
            cortest = cor.test(target_pert, source_pert)
            KP2KP$p[i] = cortest$p.value
            KP2KP$sign[i] = sign(cortest$estimate)
        }
        # using lm
        # slope = summary(lm(target_pert ~ source_pert -1))$coefficients
        # KP2KP$p[i] = slope[length(slope)]
        # KP2KP$sign[i] = sign(slope[1])
    }
}

KP2KP$p_adj = p.adjust(KP2KP$p, "BH")
KP2KP$weight = qbeta(KP2KP$p_adj, 1, 20, lower.tail=F) * KP2KP$sign

write.table(KP2KP, "KP2KP_p.tsv", sep="\t", row.names=F, quote=F)

# analysis
sum(KP2KP$p_adj < 5e-2, na.rm=T)
sum(KP2KP$p_adj < 1e-6, na.rm=T)
sum(KP2KP$p_adj < 1e-12, na.rm=T)
hist(KP2KP$weight, breaks=50)

