
# settings
options(stringsAsFactors=F)

# packages
library(harmonicmeanp)

# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
perturbation = as.matrix(read.table("../../perturbation/logFC_inner.csv", sep=",", header=T, quote=""))
TF_edges = read.table("../TF_edge_weights.tsv", sep="\t", header=T, quote="")[,1:3]
TF_edges = TF_edges[TF_edges$TF %in% colnames(perturbation),]
TFs = sort(unique(TF_edges$TF))
KPs = flatten(read.table("../KP.txt"))
KPs = KPs[KPs %in% colnames(perturbation)]
KP2TF = read.table("KP2TF_p.tsv", sep="\t", quote="", header=T)
KP2TF_sig = KP2TF[KP2TF$p < 0.05,]

KP2KP = data.frame(KP=rep(KPs, each=length(KPs)), Target=rep(KPs,length(KPs)), p=NA)
KP2KP = KP2KP[KP2KP$KP != KP2KP$Target,]  # no self-loops
for (i in seq(1,nrow(KP2KP),1)) {
    source = KP2KP$KP[i]
    target = KP2KP$Target[i]
    targetTFs = KP2TF_sig$Target[KP2TF_sig$KP == target]
    if (length(targetTFs) > 0) {
        case    = KP2TF$p[KP2TF$KP == source & (KP2TF$Target%in%targetTFs)]
        control = KP2TF$p[KP2TF$KP == target & (KP2TF$Target%in%targetTFs)]
        # https://www.statisticssolutions.com/kendalls-tau-and-spearmans-rank-correlation-coefficient/ 
        KP2KP$p[i] = cor.test(case, control, method="kendall")$p.value
    }
}

sum(KP2KP$p < 5e-2, na.rm=T)
sum(KP2KP$p < 1e-3, na.rm=T)
sum(KP2KP$p < 1e-6, na.rm=T)
sum(KP2KP$p < 1e-12, na.rm=T)
sum(KP2KP$p < 1e-15, na.rm=T)
sum(KP2KP$p < 1e-18, na.rm=T)
sum(KP2KP$p < 1e-21, na.rm=T)
hist(qbeta(na.omit(KP2KP$p), 1, 20, lower.tail=F), breaks=80)
hist(qbeta(na.omit(KP2KP$p[KP2KP$p < 0.05]), 1, 20, lower.tail=F), breaks=80)





