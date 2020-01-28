
options(stringsAsFactors=F)

# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
KP2TF = read.table("KP2TF_p.tsv", sep="\t", quote="", header=T)
KP2KP = read.table("KP2KP_p.tsv", sep="\t", quote="", header=T)
KP_edges = data.frame(KP    =c(KP2TF$KP,     KP2KP$KP),
                      Target=c(KP2TF$TF,     KP2KP$Target),
                      weight=c(KP2TF$weight, KP2KP$weight))

# adjacency matrix
KPs = flatten(read.table("../KP.txt"))
TFs = flatten(read.table("../TF.txt"))
KPTFs = c(KPs, TFs)
adjacency = as.matrix(sparseMatrix(i=match(KP_edges$Target, KPTFs), j=match(KP_edges$KP, KPs), x=KP_edges$weight, 
                                   dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
write.table(adjacency, "WP.csv", sep=",", quote=F)
write.table(adjacency, "WP.mat", sep=" ", quote=F, col.names=F, row.names=F)

adjacency_noise = adjacency
lacking = is.na(adjacency_noise) | (adjacency_noise==0)
noise_sd = 1/sqrt(sum(dim(adjacency)^2))
adjacency_noise[lacking] = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))[lacking]

write.table(adjacency_noise, "WP_noise.csv", sep=",", quote=F)
write.table(adjacency_noise, "WP_noise.mat", sep=" ", quote=F, col.names=F, row.names=F)


# analysis

KP_sig = KP_edges[abs(KP_edges$weight) > .1,]
KP_sig = paste(KP_sig$KP, KP_sig$Target)

P_eval = read.table("../../evaluation/P_eval.tsv", sep="\t", header=T, check.names=F)
gold1 = P_eval$goldstandard1 == 1 & !is.na(P_eval$goldstandard1)
gold2 = P_eval$goldstandard2 == 1 & !is.na(P_eval$goldstandard2)
gold3 = P_eval$goldstandard3 == 1 & !is.na(P_eval$goldstandard3)
gold1 = paste(P_eval$Source[gold1], P_eval$Target[gold1])
gold2 = paste(P_eval$Source[gold2], P_eval$Target[gold2])
gold3 = paste(P_eval$Source[gold3], P_eval$Target[gold3])

length(intersect(KP_sig, gold1)); length(setdiff(KP_sig, gold1)); length(setdiff(gold1, KP_sig))
length(intersect(KP_sig, gold2)); length(setdiff(KP_sig, gold2)); length(setdiff(gold2, KP_sig))
length(intersect(KP_sig, gold3)); length(setdiff(KP_sig, gold3)); length(setdiff(gold3, KP_sig))






