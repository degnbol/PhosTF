
# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/network/KP_priors")

KPs = flatten(read.table("../KP.txt"))
TF_edges = read.table("../TF_edge_weights.tsv", sep="\t", header=T, quote="")[,1:2]
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

# write.table(KP_edges, "KP2TF_all.tsv", sep="\t", quote=F, row.names=F)
# KP_edges = read.table("KP2TF_all.tsv", sep="\t", quote="", header=T)

KP_edges = KP_edges[KP_edges$p < 0.05,]
KP_edges$weight = qnorm(KP_edges$p, mean=-.35, sd=.3, lower.tail=F)
hist(KP_edges$weight, breaks=50)

write.table(KP_edges, "KP2TF.tsv", sep="\t", quote=F, row.names=F)



# adjacency matrix
KPs = flatten(read.table("../KP.txt"))
TFs = flatten(read.table("../TF.txt"))
KPTFs = c(KPs,TFs)
adjacency = as.matrix(sparseMatrix(i=match(KP_edges$Target, KPTFs), j=match(KP_edges$KP, KPs), x=KP_edges$weight, 
                                   dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
write.table(adjacency, "KP_edges.csv", sep=",", quote=F)
write.table(adjacency, "KP_edges.mat", sep=" ", quote=F, col.names=F, row.names=F)






KP_edges = paste(KP_edges$KP, KP_edges$Target)

P_eval = read.table("../../evaluation/P_eval.tsv", sep="\t", header=T, check.names=F)
gold1 = P_eval$goldstandard1 == 1 & !is.na(P_eval$goldstandard1)
gold2 = P_eval$goldstandard2 == 1 & !is.na(P_eval$goldstandard2)
gold3 = P_eval$goldstandard3 == 1 & !is.na(P_eval$goldstandard3)
gold1 = paste(P_eval$Source[gold1], P_eval$Target[gold1])
gold2 = paste(P_eval$Source[gold2], P_eval$Target[gold2])
gold3 = paste(P_eval$Source[gold3], P_eval$Target[gold3])

length(intersect(KP_edges, gold1))
length(intersect(KP_edges, gold2))
length(intersect(KP_edges, gold3))
length(setdiff(KP_edges, gold1)); length(setdiff(gold1, KP_edges))
length(setdiff(KP_edges, gold2)); length(setdiff(gold2, KP_edges))
length(setdiff(KP_edges, gold3)); length(setdiff(gold3, KP_edges))






