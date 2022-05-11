#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(here))
root = here()
d = function(relpath) {
    paste0(root, "/", relpath)
}

dt.true = fread(d("results/4-yeastNetworkReconstruction/P_edges/P_edges.tsv"))
# dt.pred = fread(d("results/4-yeastNetworkReconstruction/infer/WP_infer-inner-strict-0.0-0.1.tsv"))
dt.pred = fread(commandArgs(TRUE))


dt.true[, Ref:=NULL]
dt.true$True = TRUE

dt.pred.melt = melt(dt.pred, id.vars="_", variable.name="P", value.name="Pred", variable.factor=FALSE)
setnames(dt.pred.melt, "_", "Target")

dt = merge(dt.pred.melt, dt.true, by=c("P", "Target"), all.x=TRUE)
dt[is.na(True), True:=FALSE]

roc(dt$True, dt$Pred)$auc

