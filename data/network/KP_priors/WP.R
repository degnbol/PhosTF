
options(stringsAsFactors=F)

# functions
flatten = function(x) as.vector(as.matrix(x))

# combine pvals using fishers method https://en.wikipedia.org/wiki/Fisher%27s_method 
# chisq = -2 sum(ln(p-values)); pval = 1-pchisq(chisq, df=2length(p-values))
fisher.method.log = function(pvals) {
    df = 2*length(pvals)
    pchisq(-2*sum(log(pvals),na.rm=TRUE),df,lower.tail=FALSE,log.p=T)
}


# main

setwd("~/cwd/data/network/KP_priors")
# read
TF_edges = read.table("../TF_edge_weights.tsv", sep="\t", header=T, quote="")[,1:3]
perturbation = as.matrix(read.table("../../perturbation/logFC_inner.csv", sep=",", header=T, quote=""))
TF_edges = TF_edges[TF_edges$TF %in% colnames(perturbation),]
TFs = sort(unique(TF_edges$TF))

KP_edges = read.table("KP2TF_p.tsv", sep="\t", quote="", header=T)
KP_edges = KP_edges[KP_edges$p < 0.05,]



KP_edges$weight = qbeta(KP_edges$p, 1, 20, lower.tail=F)
hist(KP_edges$weight, breaks=50)

# find sign of KP->TF edge like this:
# either take sign of correlation of logFC values for KP and TF KOs or ...
# make sets of TF targets where TF KO and KP KO sign of target logFC agrees and set where it doesn't agree.
# combine p-vals for agree vs. disagree sets and decide sign of KP->TF edge from which is most significant. 

KP_edges$corsign = NA
KP_edges$fishsign = NA
for (i in 1:nrow(KP_edges)) {
    KP = KP_edges$KP[i]
    TF = KP_edges$Target[i]
    TF_targets = TF_edges[TF_edges$TF == TF,]
    KP_edges$corsign[i] = sign(cor(perturbation[TF_targets$Target,KP], perturbation[TF_targets$Target,TF]))
    agree = sign(perturbation[TF_targets$Target,KP]) == sign(perturbation[TF_targets$Target,TF])
    p_agree = fisher.method.log(TF_targets$Pval[agree])
    p_disagree = fisher.method.log(TF_targets$Pval[!agree])
    if (p_agree < p_disagree) {KP_edges$fishsign[i] = +1}
    if (p_agree > p_disagree) {KP_edges$fishsign[i] = -1}
}

sum(is.na(KP_edges$corsign))
sum(is.na(KP_edges$fishsign))
sum(is.na(KP_edges$corsign) & is.na(KP_edges$fishsign))
sum(KP_edges$corsign != KP_edges$fishsign, na.rm=T) / nrow(na.omit(KP_edges))
# ohno


write.table(KP_edges, "KP2TF.tsv", sep="\t", quote=F, row.names=F)


# adjacency matrix
KPs = flatten(read.table("../KP.txt"))
TFs = flatten(read.table("../TF.txt"))
KPTFs = c(KPs,TFs)
adjacency = as.matrix(sparseMatrix(i=match(KP_edges$Target, KPTFs), j=match(KP_edges$KP, KPs), x=KP_edges$weight, 
                                   dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
write.table(adjacency, "WP.csv", sep=",", quote=F)
write.table(adjacency, "WP.mat", sep=" ", quote=F, col.names=F, row.names=F)







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






