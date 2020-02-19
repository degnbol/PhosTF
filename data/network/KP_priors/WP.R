
library(Matrix)
library(fdrtool)

options(stringsAsFactors=F)

# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
KPs = flatten(read.table("../KP.txt"))
TFs = flatten(read.table("../TF.txt"))
KPTFs = c(KPs, TFs)
KP2TF = read.table("KP2TF.tsv", sep="\t", quote="", header=T)
KP2KP = read.table("archive/KP2KP_p.tsv", sep="\t", quote="", header=T)
colnames(KP2TF)[colnames(KP2TF)=="TF"] = "Target"
KP2KP = KP2KP[,colnames(KP2KP)!="p_adj"]
KP_edges = rbind(KP2TF, KP2KP)

KP2TF$q = fdrtool(KP2TF$p, statistic="pvalue", plot=FALSE)$qval
KP2TF$normlogp = log10(KP2TF$p) / min(log10(KP2TF$p)) * KP2TF$sign
KP2TF_FDR10 = KP2TF[KP2TF$q < .1,]
KP2TF_FDR20 = KP2TF[KP2TF$q < .2,]

KP2TF$gauss = qhalfnorm(KP2TF$p, lower.tail=FALSE)
KP2TF$gauss01 = qhalfnorm(KP2TF$p, theta=sd2theta(.1), lower.tail=FALSE)
KP2TF$gauss[KP2TF$gauss == Inf] = max(KP2TF$gauss[KP2TF$gauss != Inf])
KP2TF$gauss01[KP2TF$gauss01 == Inf] = max(KP2TF$gauss01[KP2TF$gauss01 != Inf])
KP2TF$gauss = KP2TF$gauss*KP2TF$sign
KP2TF$gauss01 = KP2TF$gauss01*KP2TF$sign

# adjacency matrix
adjacency = as.matrix(sparseMatrix(i=match(KP_edges$Target, KPTFs), j=match(KP_edges$KP, KPs), x=KP_edges$weight, 
                                   dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF = as.matrix(sparseMatrix(i=match(KP2TF$Target, KPTFs), j=match(KP2TF$KP, KPs), x=KP2TF$weight, 
                                   dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF_FDR10 = as.matrix(sparseMatrix(i=match(KP2TF_FDR10$Target, KPTFs), j=match(KP2TF_FDR10$KP, KPs), x=KP2TF_FDR10$sign, 
                                               dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF_FDR20 = as.matrix(sparseMatrix(i=match(KP2TF_FDR20$Target, KPTFs), j=match(KP2TF_FDR20$KP, KPs), x=KP2TF_FDR20$sign, 
                                               dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF_log = as.matrix(sparseMatrix(i=match(KP2TF$Target, KPTFs), j=match(KP2TF$KP, KPs), x=KP2TF$normlogp, 
                                         dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF_gauss = as.matrix(sparseMatrix(i=match(KP2TF$Target, KPTFs), j=match(KP2TF$KP, KPs), x=KP2TF$gauss, 
                                             dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))
adjacency_KP2TF_gauss01 = as.matrix(sparseMatrix(i=match(KP2TF$Target, KPTFs), j=match(KP2TF$KP, KPs), x=KP2TF$gauss01, 
                                               dims=list(length(KPTFs), length(KPs)), dimnames=list(KPTFs, KPs)))


adjacency[is.na(adjacency)] = 0

get_adjacency_noise = function(adjacency) {
    adjacency_noise = adjacency
    noise_sd = 1/sqrt(prod(dim(adjacency)))
    lacking = adjacency_noise==0
    adjacency_noise[lacking] = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))[lacking]
    adjacency_noise
}



write.table(adjacency, "WP.csv", sep=",", quote=F)
write.table(adjacency, "WP.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacency_KP2TF, "WP_KP2TF.csv", sep=",", quote=F)
write.table(adjacency_KP2TF, "WP_KP2TF.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency), "WP_noise.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF), "WP_noise_KP2TF.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(sign(adjacency), "WP_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(sign(adjacency_KP2TF), "WP_KP2TF_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(sign(adjacency)), "WP_noise_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(sign(adjacency_KP2TF)), "WP_noise_KP2TF_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacency_KP2TF_FDR10, "WP_FDR10_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacency_KP2TF_FDR20, "WP_FDR20_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF_FDR10), "WP_noise_FDR10_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF_FDR20), "WP_noise_FDR20_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacency_KP2TF_log, "WP_KP2TF_log.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF_log), "WP_noise_KP2TF_log.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF_gauss), "WP_noise_KP2TF_gauss.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(get_adjacency_noise(adjacency_KP2TF_gauss01), "WP_noise_KP2TF_gauss01.mat", sep=" ", quote=F, col.names=F, row.names=F)

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






