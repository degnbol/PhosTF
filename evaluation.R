#!/usr/bin/env Rscript

# packages
library(reshape2)
library(Matrix)
suppressPackageStartupMessages(library(pROC))

# functions
flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
read.matrix = function(x) as.matrix(read.table(x))
melt_matrix = function(x) {
    out = melt(as.matrix(x))
    colnames(out) = c("rownames", "colnames", "value")
    out
}


setwd("/Users/christian/GoogleDrev/PKTFX")

# read
P_fname = "data/evaluation/P_eval.tsv"
T_fname = "data/evaluation/T_eval.tsv"
WP_fname = "data/inference/05/WP_infer_1.mat"
WT_fname = "data/inference/05/WT_infer_1.mat"
WT_prior_fname = "data/pertubation/WT_prior.mat"
KP_fname = "data/pertubation/KP.txt"
TF_fname = "data/pertubation/TF.txt"
PTX_fname = "data/pertubation/PTX.txt"
P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
T_eval = read.table(T_fname, header=T, sep="\t", quote="", check.names=F)
WP = read.matrix(WP_fname)
WT = read.matrix(WT_fname)
WT_prior = read.matrix(WT_prior_fname)
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
PTX = read.vector(PTX_fname)
PT = c(KP,TF)
colnames(WP) = KP
colnames(WT) = TF
colnames(WT_prior) = TF
rownames(WP) = c(KP, TF)
rownames(WT) = PTX
rownames(WT_prior) = PTX
# swap rownames and colnames columns so the order will be source then target
P_edges = melt_matrix(WP)[,c(2,1,3)]
T_edges = melt_matrix(WT)[,c(2,1,3)]
T_prior_edges = melt_matrix(WT_prior)[,c(2,1,3)]
colnames(P_edges) = c("P",  "Target", "marker")
colnames(T_edges) = c("TF", "Target", "marker")
colnames(T_prior_edges) = c("TF", "Target", "mask")
# sanity check
all(P_edges$P == P_eval$Source)
all(P_edges$Target == P_eval$Target)
all(T_edges$TF == T_eval$Source)
all(T_edges$Target == T_eval$Target)
all(T_prior_edges$TF == T_eval$Source)
all(T_prior_edges$Target == T_eval$Target)

# mask with prior
T_edges_masked = T_edges[T_prior_edges$mask == 1,]
T_eval_masked = T_eval[T_prior_edges$mask == 1,]

# hold all data in eval sets
P_eval$marker = P_edges$marker
T_eval$marker = T_edges$marker
T_eval_masked$marker = T_edges_masked$marker

# remove diagonals
P_eval = P_eval[as.character(P_eval$Source) != as.character(P_eval$Target),]
T_eval = T_eval[as.character(T_eval$Source) != as.character(T_eval$Target),]
T_eval_masked = T_eval_masked[as.character(T_eval_masked$Source) != as.character(T_eval_masked$Target),]

cor(P_eval$yeastkid[!is.na(P_eval$yeastkid)], abs(P_eval$marker[!is.na(P_eval$yeastkid)]))
cor(P_eval$reaction[!is.na(P_eval$reaction)], abs(P_eval$marker[!is.na(P_eval$reaction)]))
cor(P_eval$ptmod[!is.na(P_eval$ptmod)], abs(P_eval$marker[!is.na(P_eval$ptmod)]))
cor(P_eval$inhibition[!is.na(P_eval$inhibition)], -(P_eval$marker[!is.na(P_eval$inhibition)]))
cor(P_eval$activation[!is.na(P_eval$activation)], +(P_eval$marker[!is.na(P_eval$activation)]))
cor(P_eval$expression[!is.na(P_eval$expression)], abs(P_eval$marker[!is.na(P_eval$expression)]))
cor(P_eval$catalysis[!is.na(P_eval$catalysis)], abs(P_eval$marker[!is.na(P_eval$catalysis)]))
cor(P_eval$n_datasets, abs(P_eval$marker))
auc(roc(P_eval$n_datasets > 0, abs(P_eval$marker)))
auc(roc(P_eval$parca, abs(P_eval$marker)))
auc(roc(P_eval$fasolo, abs(P_eval$marker)))
auc(roc(P_eval$biogrid != "", abs(P_eval$marker)))
auc(roc(P_eval$biogrid == "kinase", +(P_eval$marker)))
auc(roc(P_eval$biogrid == "phosphatase", -(P_eval$marker)))


evaluate_T = function(dataset) {
    out = c()
    cor_names = c("n_datasets", "reaction", "expression", "catalysis", "workman SLL", "workman pval")
    cor_names_pos = c("activation")
    cor_names_neg = c()
    for (name in cor_names) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], abs(dataset$marker[valid]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    for (name in cor_names_pos) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], +(dataset$marker[valid]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    for (name in cor_names_neg) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], -(dataset$marker[valid]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    paste(out, collapse="\n")
}

evaluate_auc_T = function(dataset) {
    aucs = c(
        auc(roc(dataset$`yeastract expression` != "", abs(dataset$marker))),
        auc(roc(dataset$`yeastract expression` == "activator", +(dataset$marker))),
        auc(roc(dataset$`yeastract expression` == "inhibitor", -(dataset$marker))),
        auc(roc(!is.na(dataset$`yeastract binding`), abs(dataset$marker))),
        auc(roc(!is.na(dataset$balaji), abs(dataset$marker)))
    )
    paste(c("yeastract expression", "yeastract activator", "yeastract inhibitor", "yeastract binding", "balaji"), aucs, sep="\t", collapse="\n")
}


cat(evaluate_T(T_eval_masked), "\n")
aucs_T = evaluate_auc_T(T_eval_masked)
cat(aucs_T, "\n")

# stratification
# script



