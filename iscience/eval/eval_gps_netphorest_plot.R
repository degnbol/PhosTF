#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
library(ggplot2)
library(pammtools)


setwd("~/PhosTF/iscience/eval/")
TPRs.FPRs = fread("roc.tsv")


ggplot() + 
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method))




ROC.dt = data.table()
for (meth in unique(TPRs.FPRs$method)) {
    ROC.dt = rbind(ROC.dt, merge(
        TPRs.FPRs[(method==meth) & (beta == "beta_est")],
        CJ(FPR=TPRs.FPRs[(method==meth) & (beta == "beta_est"), FPR], bound=c("upper", "lower"), method=meth, beta="beta_est", unique=T)
        , by=c("FPR", "bound", "method", "beta"), all=T),
        merge(
            TPRs.FPRs[(method==meth) & (beta != "beta_est")],
            CJ(FPR=TPRs.FPRs[(method==meth) & (beta != "beta_est"), FPR], bound=c("upper", "lower"), method=meth, beta=c("beta_upper", "beta_lower"), unique=T)
            , by=c("FPR", "bound", "method", "beta"), all=T)
    )
}

# given a vector values with some NA values, copy the last value that is not NA, replacing NA
fill_last_nona = function(vals) {
    while(any(is.na(vals))) {
        last = c(vals[1], vals[-length(vals)])
        vals[is.na(vals)] = last[is.na(vals)]
    }
    vals
}

ROC.dt[, TPR:=fill_last_nona(TPR), by=c("bound", "method", "beta")]
ROC.dt = ROC.dt[order(FPR)]
ROC.dt.bounds = merge(
    unique(ROC.dt[bound=="upper", .(TPR_ub=max(TPR)), by=c("FPR","method","beta")]),
    unique(ROC.dt[bound=="lower", .(TPR_lb=max(TPR)), by=c("FPR","method","beta")]), 
    by=c("FPR","method","beta"))


ggplot() +
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_stepribbon(data=ROC.dt.bounds[beta=="beta_est"], mapping=aes(x=FPR, ymin=TPR_lb, ymax=TPR_ub, fill=method), alpha=0.4)

ggplot() +
    geom_step(data=TPRs.FPRs[(beta=="beta_lower") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_lower") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_stepribbon(data=ROC.dt.bounds[beta=="beta_lower"], mapping=aes(x=FPR, ymin=TPR_lb, ymax=TPR_ub, fill=method), alpha=0.3) +
    geom_step(data=TPRs.FPRs[(beta=="beta_upper") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_upper") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_stepribbon(data=ROC.dt.bounds[beta=="beta_upper"], mapping=aes(x=FPR, ymin=TPR_lb, ymax=TPR_ub, fill=method), alpha=0.3)


ROC.dt.lim = ROC.dt.bounds[beta!="beta_est", .(TPR_ub=max(TPR_ub), TPR_lb=min(TPR_lb)), by=c("method", "FPR")]

ggplot() +
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_est") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_stepribbon(data=ROC.dt.lim, mapping=aes(x=FPR, ymin=TPR_lb, ymax=TPR_ub, fill=method), alpha=0.25)





