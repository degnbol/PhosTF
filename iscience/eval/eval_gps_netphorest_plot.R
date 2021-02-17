#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
library(ggplot2)
library(pammtools)


setwd("~/PhosTF/iscience/eval/")
TPRs.FPRs = fread("roc.tsv")

TPRs.FPRs[(method=="PhosTF") & (bound=="median") & (beta=="beta_lower") & (P >= 8776)]
TPRs.FPRs[(method=="PhosTF") & (bound=="median") & (beta=="beta_upper") & (P >= 8776)]


ggplot() + 
    geom_step(data=TPRs.FPRs[(beta=="beta_lower") & (bound=="upper")], mapping=aes(x=FPR, y=TPR, color=method)) +
    geom_step(data=TPRs.FPRs[(beta=="beta_lower") & (bound=="lower")], mapping=aes(x=FPR, y=TPR, color=method))




ROC.dt = data.table()
for (meth in unique(TPRs.FPRs$method)) {
    ROC.dt = rbind(ROC.dt, merge(TPRs.FPRs[method==meth],
            CJ(FPR=TPRs.FPRs[method==meth, FPR], bound=c("upper", "lower"), method=meth, beta=c("beta_upper", "beta_lower"), unique=T)
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


ribbon_alpha = 0.4
ribbon_stroke_alpha = 0.8

ggplot() + 
    geom_path(data=data.table(x=c(0,1), y=c(0,1)), mapping=aes(x=x, y=y), alpha=0.5) +
    theme_linedraw() + scale_x_continuous(name="FPR", expand=c(0,0)) + scale_y_continuous(name="TPR", expand=c(0,0)) + 
    theme(panel.grid.minor=element_line(color="lightgray"), panel.spacing = unit(2, "lines")) + # panel spacing separates the two subplots a bit more
    facet_grid(cols=vars(beta)) +
    geom_stepribbon(data=ROC.dt.bounds, mapping=aes(x=FPR, ymin=TPR_lb, ymax=TPR_ub, fill=method), alpha=ribbon_alpha) +
    geom_step(data=TPRs.FPRs[(bound=="median") & (method=="PhosTF")], mapping=aes(x=FPR, y=TPR), color="blue") +
    coord_fixed()  # force square


ggsave("roc_boot.pdf", width=11, height=5)




