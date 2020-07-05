#!/usr/bin/env Rscript
library(data.table)
library(Matrix)
library(fdrtool)

options(stringsAsFactors=F)

# functions
flatten = function(x) as.vector(as.matrix(x))

# main

setwd("~/cwd/data/network/KP_priors")
# read
KPs = fread("../KP_protein.tsv")$ORF
TFs = fread("../TF.txt", header=F)$V1
KPTFs = c(KPs, TFs)
#KP2KP = fread("KP2KP.tsv", sep="\t", quote="", header=T)
KP2TF = fread("KP2TF.tsv", sep="\t", quote="", header=T)
colnames(KP2TF)[colnames(KP2TF)=="TF"] = "substrate"
nKP = length(KPs)
nTF = length(TFs)
nPT = nKP+nTF

MP = fread("../KP_mask.mat", sep=" ", header=F)
noise_sd = sqrt(1/(sum(MP)-length(KPs))) # diagonal should be zero

KP2TF_FDR10 = KP2TF[q < .1,]
KP2TF_FDR20 = KP2TF[q < .2,]
#KP2KP_FDR20 = KP2KP[q < .2,]


#KP2KPTF_FDR20 = rbind(KP2KP_FDR20[,c("KP", "substrate", "sign", "median_weight")],
#                      KP2TF_FDR20[,c("KP", "substrate", "sign", "median_weight")])


setnames(KP2TF, "p.value", "p")
KP2TF[, normlogp:=log10(p)/min(log10(p))*sign]

KP2TF$gauss = qhalfnorm(KP2TF$p, lower.tail=FALSE)
KP2TF$gauss01 = qhalfnorm(KP2TF$p, theta=sd2theta(.1), lower.tail=FALSE)
KP2TF$gauss001 = qhalfnorm(KP2TF$p, theta=sd2theta(.01), lower.tail=FALSE)
KP2TF$gauss003 = qhalfnorm(KP2TF$p, theta=sd2theta(noise_sd), lower.tail=FALSE)
KP2TF$gauss[KP2TF$gauss == Inf] = max(KP2TF$gauss[KP2TF$gauss != Inf])
KP2TF$gauss01[KP2TF$gauss01 == Inf] = max(KP2TF$gauss01[KP2TF$gauss01 != Inf])
KP2TF$gauss001[KP2TF$gauss001 == Inf] = max(KP2TF$gauss001[KP2TF$gauss001 != Inf])
KP2TF$gauss003[KP2TF$gauss003 == Inf] = max(KP2TF$gauss003[KP2TF$gauss003 != Inf])
KP2TF$gauss = KP2TF$gauss*KP2TF$sign
KP2TF$gauss01 = KP2TF$gauss01*KP2TF$sign
KP2TF$gauss001 = KP2TF$gauss001*KP2TF$sign
KP2TF$gauss003 = KP2TF$gauss003*KP2TF$sign

# adjacency matrices
sparsematrix = function(i, j, x) {
    i = match(i, KPTFs)
    j = match(j, KPs)
    dims_ = list(length(KPTFs), length(KPs))
    dimnames_ = list(KPTFs, KPs)
    as.matrix(sparseMatrix(i=i, j=j, x=x, dims=dims_, dimnames=dimnames_))
}

add_noise = function(adjacency, noise_sd = 1/sqrt(prod(dim(adjacency)))) {
    cat(noise_sd, "\n")
    adjacency_noise = adjacency
    lacking = adjacency_noise==0
    adjacency_noise[lacking] = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))[lacking]
    adjacency_noise
}

#adjacency_median = sparsematrix(KP2KPTF$substrate, KP2KPTF$KP, KP2KPTF$median_weight)
#adjacency_sign = sparsematrix(KP2KPTF$substrate, KP2KPTF$KP, KP2KPTF$sign)

adjacency_KP2TF = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$median_weight)
adjacency_KP2TF_FDR10_median = sparsematrix(KP2TF_FDR10$substrate, KP2TF_FDR10$KP, KP2TF_FDR10$median_weight) 
adjacency_KP2TF_FDR20_median = sparsematrix(KP2TF_FDR20$substrate, KP2TF_FDR20$KP, KP2TF_FDR20$median_weight) 
adjacency_KP2TF_FDR10_sign = sparsematrix(KP2TF_FDR10$substrate, KP2TF_FDR10$KP, KP2TF_FDR10$sign) 
adjacency_KP2TF_FDR20_sign = sparsematrix(KP2TF_FDR20$substrate, KP2TF_FDR20$KP, KP2TF_FDR20$sign) 
adjacency_KP2TF_log = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$normlogp)
adjacency_KP2TF_gauss = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$gauss)
adjacency_KP2TF_gauss01 = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$gauss01)
adjacency_KP2TF_gauss001 = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$gauss001)
adjacency_KP2TF_gauss003 = sparsematrix(KP2TF$substrate, KP2TF$KP, KP2TF$gauss003)



#fwrite(adjacency_median, "WP_median_FDR20.mat", sep=" ", row.names=F, col.names=F)
#fwrite(add_noise(adjacency_median), "WP_median_FDR20_noise.mat", sep=" ", row.names=F, col.names=F)
#fwrite(adjacency_sign, "WP_sign_FDR20.mat", sep=" ", row.names=F, col.names=F)
#fwrite(add_noise(adjacency_sign), "WP_sign_FDR20_noise.mat", sep=" ", row.names=F, col.names=F)


fwrite(adjacency_KP2TF_FDR20_median, "WP_median_KP2TF_FDR20.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(adjacency_KP2TF_FDR20_median), "WP_median_KP2TF_FDR20_noise.mat", sep=" ", row.names=F, col.names=F)
fwrite(sign(adjacency_KP2TF_FDR20_median), "WP_sign_KP2TF_FDR20.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(sign(adjacency_KP2TF_FDR20_median)), "WP_sign_KP2TF_FDR20_noise.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(adjacency_KP2TF_gauss, 1), "WP_KP2TF_gauss.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(adjacency_KP2TF_gauss01, .1), "WP_KP2TF_gauss01.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(adjacency_KP2TF_gauss001, .01), "WP_KP2TF_gauss001.mat", sep=" ", row.names=F, col.names=F)
fwrite(add_noise(adjacency_KP2TF_gauss003, noise_sd), "WP_KP2TF_gauss003.mat", sep=" ", row.names=F, col.names=F)


fwrite(matrix(rnorm(nPT*nKP, sd=0.01), nrow=nPT, ncol=nKP),
       "WP_gauss001.mat", sep=" ", row.names=F, col.names=F)



