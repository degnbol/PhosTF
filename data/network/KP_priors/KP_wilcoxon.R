
# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
KPs = flatten(read.table("../KP.txt"))
TF_edges = read.table("../TF_edge_weights.tsv", sep="\t", header=T, quote="")[,1:3]
perturbation = as.matrix(read.table("../../perturbation/logFC_inner.csv", sep=",", header=T, quote=""))
KPs = KPs[KPs %in% colnames(perturbation)]
TF_edges = TF_edges[TF_edges$TF %in% colnames(perturbation),]
TFs = sort(unique(TF_edges$TF))

KP_edges = data.frame(KP=rep(KPs, each=length(TFs)), Target=rep(TFs,length(KPs)), p=NA)
for (i in 1:nrow(KP_edges)) {
    KP = KP_edges$KP[i]
    TF = KP_edges$Target[i]
    TF_targets = TF_edges$Target[TF_edges$TF == TF]
    case = abs(perturbation[TF_targets,KP])
    control = abs(perturbation[!(rownames(perturbation) %in% TF_targets),KP])
    KP_edges$p[i] = wilcox.test(case, control, "greater")$p.value
}

write.table(KP_edges, "KP2TF_all.tsv", sep="\t", quote=F, row.names=F)


