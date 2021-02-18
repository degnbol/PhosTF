#!/usr/bin/env Rscript
# eval gps, netphorest and phosTF on evaluation data with precision, recall, specificity, auc, roc, etc.
suppressPackageStartupMessages(library(data.table))
library(matrixStats)

setwd("~/PhosTF/iscience/eval/")

DT = fread("../../data/evaluation/P_eval.tsv")
DT[,biogrid:=biogrid!=""]
DT[,fiedler:=fiedler!=""]
DT[,fasolo:=!is.na(fasolo)]
DT[,parca:=!is.na(parca)]
DT[,ptacek:=!is.na(ptacek)]
DT[,yeastkid05:=yeastkid>4.52 & !is.na(yeastkid)]
DT[,yeastkid01:=yeastkid>6.40 & !is.na(yeastkid)]
DT[,ptmod5:=ptmod>500 & !is.na(ptmod)]
DT[,eval:=biogrid|fiedler|fasolo|parca|ptacek|yeastkid05|ptmod5]
setnames(DT, c("Source", "Target"), c("kinase", "substrate"))
GPS = fread("../gps/gps_edge_scores.tsv")
setnames(GPS, "score_max_norm", "gps")
DT = merge(DT, GPS, all=T, by=c("kinase", "substrate"))
PhosTF = fread("../../data/inference/KP_edges.tsv", drop="q")
setnames(PhosTF, "marker", "phostf")
DT = merge(DT, PhosTF, all=T, by=c("kinase", "substrate"))

# discard everything we are not going to use for better overview
DT = DT[!is.na(eval), .(eval, phostf, netphorest, gps)]
# split for each method
pho = DT[!is.na(phostf), .(eval, score=phostf, rank=rank(-phostf))]
net = DT[!is.na(netphorest), .(eval, score=netphorest, rank=rank(-netphorest))]
gps = DT[!is.na(gps), .(eval, score=gps, rank=rank(-gps))]

# all unlabeled are considered negatives here, i.e. beta = 0
# pho[,plot(roc(eval, score, direction=">"))]
# net[,plot(roc(eval, score, direction="<"))]
# gps[,plot(roc(eval, score, direction="<"))]


## estimate CI with bootstrap from the known positives

get_boot_tpr_quantiles = function(P.L.ranks, rank.range, n.replicates=1, probs=c(0.025, 0.975)) {
    P.L.ranks.boot = matrix(sample(P.L.ranks, size=length(P.L.ranks) * n.replicates, replace=T), ncol=n.replicates)
    sapply(rank.range, function (r) quantile(colMeans(P.L.ranks.boot <= r), probs=probs))
}

get_T = function(P.L.ranks, U.ranks, n.replicates=100) {
    rank.range = sort(unique(c(P.L.ranks, U.ranks)))
    
    TPR.L.CI = get_boot_tpr_quantiles(P.L.ranks, rank.range, n.replicates)
    
    # look at bootstrapped CI of rank vs TPR
    # plot(rank.range, rowMeans(tprs.boot), type='l')
    # lines(rank.range, TPR.L.CI[1,])
    # lines(rank.range, TPR.L.CI[2,])
    list(lb=TPR.L.CI[1,], ub=TPR.L.CI[2,])
}

# eq. 10
get_theta.ub = function(T.ub, n.P.U.star) {
    ceiling(T.ub * n.P.U.star)
}
# eq. 13
get_theta.lb = function(T.lb, n.P.U.star) {
    floor(T.lb * n.P.U.star)
}


get_n.head = function(ranks, rank.range) {
    rowSums(matrix(ranks, nrow=length(rank.range), ncol=length(ranks), byrow=T) <= matrix(rank.range, nrow=length(rank.range), ncol=length(ranks), byrow=F))
}
get_n.tail = function(ranks, rank.range) {
    rowSums(matrix(ranks, nrow=length(rank.range), ncol=length(ranks), byrow=T) >  matrix(rank.range, nrow=length(rank.range), ncol=length(ranks), byrow=F))
}

# eq. 11
get_n.head_P.U.star = function(n.head_U, n.tail_U, rank.range, n.P.U.star, theta) {
    cond = n.P.U.star - theta <= n.tail_U
    rowMins(cbind(n.head_U, theta)) * cond + (n.P.U.star - n.tail_U) * (1-cond)
}


get_contingency = function(P.L.ranks, U.ranks, n.P.U.star, theta) {
    rank.range = sort(unique(c(P.L.ranks, U.ranks)))
    n.U = length(U.ranks)
    n.P.L = length(P.L.ranks)
    n.N.U.star = n.U - n.P.U.star
    n.head_U = get_n.head(U.ranks, rank.range)
    n.tail_U = n.U - n.head_U
    
    # these are the 4 values for the contingency table in eq. 12 for any rank r:
    n.head_P.U.star = get_n.head_P.U.star(n.head_U, n.tail_U, rank.range, n.P.U.star, theta)
    n.tail_P.U.star = n.P.U.star - n.head_P.U.star
    n.head_N.U.star = n.head_U - n.head_P.U.star
    n.tail_N.U.star = n.N.U.star - n.head_N.U.star
    
    # the 4 values for the contingency table of known labels
    n.head_P.L = get_n.head(P.L.ranks, rank.range)
    n.tail_P.L = n.P.L - n.head_P.L
    n.head_N.L = 0  # no known negatives
    n.tail_N.L = 0  # no known negatives
    
    # contingency table from start of sec. 4.1
    list(
        TP = n.head_P.L + n.head_P.U.star,
        FN = n.tail_P.L + n.tail_P.U.star,
        FP = n.head_N.L + n.head_N.U.star,
        TN = n.tail_N.L + n.tail_N.U.star
    )
}


get_TPRs.FPRs = function(P.L.ranks, U.ranks, beta, n.replicates=2000) {
    n.P.L = length(P.L.ranks)
    n.U = length(U.ranks)
    # not sure if need to round. It makes sense to be a whole number but the result seems to look the same if we don't round.
    n.P.U.star = round(n.U * beta)
    n.N.U.star = n.U - n.P.U.star
    # correct the input beta
    beta = n.P.U.star / n.U
    stopifnot(!is.na(beta))
    
    T. = get_T(P.L.ranks, U.ranks, n.replicates)
    
    con_UB = get_contingency(P.L.ranks, U.ranks, n.P.U.star, get_theta.ub(T.$ub, n.P.U.star))
    con_LB = get_contingency(P.L.ranks, U.ranks, n.P.U.star, get_theta.lb(T.$lb, n.P.U.star))
    con_med = get_contingency(P.L.ranks, U.ranks, n.P.U.star, T.$med * n.P.U.star)
    
    UB = data.table(
        TPR=con_UB$TP / (n.P.L + n.P.U.star),
        FPR=con_UB$FP / (0     + n.N.U.star),  # no known negatives
        P=con_UB$TP+con_UB$FP  # also keep the number of edges classified as a positive at a given threshold
    )
    UB$bound = "upper"
    LB = data.table(
        TPR=con_LB$TP / (n.P.L + n.P.U.star),
        FPR=con_LB$FP / (0     + n.N.U.star),   # no known negatives
        P=con_LB$TP+con_LB$FP  # also keep the number of edges classified as a positive at a given threshold
    )
    LB$bound = "lower"
    MED = data.table(
        TPR=con_med$TP / (n.P.L + n.P.U.star),
        FPR=con_med$FP / (0     + n.N.U.star),   # no known negatives
        P=con_med$TP+con_med$FP  # also keep the number of edges classified as a positive at a given threshold
    )
    MED$bound = "median"
    rbind(UB, LB, MED)
}

# min and max number of TFs that could have mediated KP->TF effects based off of 
n.known = DT[, sum(eval)]
beta_lower = (3067 - n.known) / (163 * 213)
beta_upper = (14026 - n.known) / (163 * 213)

n.replicates = 100

TPRs.FPRs = data.table()
for(beta in c(beta_lower, beta_upper)) {
    TPRs.FPRs_pho = get_TPRs.FPRs(pho[eval==T, rank], pho[eval==F, rank], beta, n.replicates)
    TPRs.FPRs_net = get_TPRs.FPRs(net[eval==T, rank], net[eval==F, rank], beta, n.replicates)
    TPRs.FPRs_gps = get_TPRs.FPRs(gps[eval==T, rank], gps[eval==F, rank], beta, n.replicates)
    TPRs.FPRs_pho$method = "PhosTF"
    TPRs.FPRs_net$method = "NetPhorest"
    TPRs.FPRs_gps$method = "GPS"
    TPRs.FPRs_beta = rbind(TPRs.FPRs_pho, TPRs.FPRs_net, TPRs.FPRs_gps)
    TPRs.FPRs_beta$beta = beta
    TPRs.FPRs = rbind(TPRs.FPRs, TPRs.FPRs_beta)
}


# expand so we have multiple TPRs for the same FPRs for comparisons that needs to be done.
TPRs.FPRs = unique(rbind(TPRs.FPRs[, .(FPR=0, TPR=0, P=0), by=c("method", "bound", "beta")], TPRs.FPRs))[order(FPR,TPR,P)]


TPRs.FPRs[beta==beta_upper, Beta:="beta_upper"]
TPRs.FPRs[beta==beta_lower, Beta:="beta_lower"]
TPRs.FPRs[, beta:=NULL]
setnames(TPRs.FPRs, "Beta", "beta")

fwrite(TPRs.FPRs, "roc.tsv", sep='\t')
