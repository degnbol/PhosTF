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


example = function() {
    
    setwd("/Users/christian/GoogleDrev/PKTFX")
    
    # read
    WP_fname = "data/inference/05/WP_infer_1.mat"
    WT_fname = "data/inference/05/WT_infer_1.mat"
    WT_prior_fname = "data/perturbation/WT_prior.mat"
    P_fname = "data/evaluation/P_eval.tsv"
    T_fname = "data/evaluation/T_eval.tsv"
    KP_fname = "data/perturbation/KP.txt"
    TF_fname = "data/perturbation/TF.txt"
    PTX_fname = "data/perturbation/PTX.txt"
    
}

# assuming when run as a script that we are in a folder such as data/inference/??/
# read
args = commandArgs(trailingOnly=T)
WP_fname = args[1]
WT_fname = args[2]
P_fname = "../../evaluation/P_eval.tsv"
T_fname = "../../evaluation/T_eval.tsv"
WT_prior_fname = "../../perturbation/WT_prior.mat"
KP_fname = "../../perturbation/KP.txt"
TF_fname = "../../perturbation/TF.txt"
PTX_fname = "../../perturbation/PTX.txt"
# end


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


# stratification
# summary(P_eval$undirected, na.rm=T)
# summary(P_eval$yeastkid, na.rm=T)
# summary(P_eval$ptmod, na.rm=T)
# summary(P_eval$expression, na.rm=T)

# sum(P_eval$yeastkid > 4.52, na.rm=T)
# sum(P_eval$yeastkid < 0, na.rm=T)
# sum(P_eval$undirected > 950, na.rm=T)
# sum(P_eval$undirected < 155, na.rm=T)
# sum(P_eval$ptmod > 300, na.rm=T)
# sum(P_eval$ptmod < 300, na.rm=T)

positives1 = (P_eval$yeastkid > 4.52) | (P_eval$undirected > 950) | (P_eval$ptmod > 300) | P_eval$fasolo
negatives1 = (P_eval$yeastkid < 0.00) | (P_eval$undirected < 155) | (P_eval$ptmod < 300)
positives2 = (P_eval$yeastkid > 4.52) | (P_eval$ptmod > 300)
negatives2 = (P_eval$yeastkid < 0.00) | (P_eval$ptmod < 300)
positives3 = P_eval$ptmod > 300
negatives3 = P_eval$ptmod < 300
P_eval$goldstandard1 = NA
P_eval$goldstandard2 = NA
P_eval$goldstandard3 = NA
P_eval$goldstandard1[positives1] = 1
P_eval$goldstandard1[negatives1] = 0
P_eval$goldstandard2[positives2] = 1
P_eval$goldstandard2[negatives2] = 0
P_eval$goldstandard3[positives3] = 1
P_eval$goldstandard3[negatives3] = 0




evaluate_cor = function(dataset, cor_names, cor_names_pos, cor_names_neg) {
    out = c()
    for (name in cor_names) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], abs(dataset[valid, "marker"]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    for (name in cor_names_pos) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], +(dataset[valid, "marker"]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    for (name in cor_names_neg) {
        valid = !is.na(dataset[name])
        PCC = cor(dataset[valid,name], -(dataset[valid, "marker"]))
        out = c(out, sprintf(paste(name, "cor:\t%.4f"), PCC))
    }
    paste(out, collapse="\n")
}


evaluate_P = function(dataset) {
    cor_names = c("yeastkid", "reaction", "ptmod", "expression", "catalysis", "netphorest", "networkin", "networkin_biogrid", "undirected", "n_datasets")
    cor_names_pos = c("activation")
    cor_names_neg = c("inhibition")
    evaluate_cor(dataset, cor_names, cor_names_pos, cor_names_neg)
}


evaluate_T = function(dataset) {
    cor_names = c("n_datasets", "reaction", "expression", "catalysis", "workman SLL", "workman pval", "undirected")
    cor_names_pos = c("activation")
    cor_names_neg = c()
    evaluate_cor(dataset, cor_names, cor_names_pos, cor_names_neg)
}


evaluate_auc_P = function(dataset) {
    aucs = c(
        auc(roc(dataset$n_datasets > 0, abs(dataset$marker))),
        auc(roc(!is.na(dataset$parca), abs(dataset$marker))),
        auc(roc(!is.na(dataset$fasolo), abs(dataset$marker))),
        auc(roc(dataset$biogrid != "", abs(dataset$marker))),
        auc(roc(dataset$biogrid == "kinase", +(dataset$marker))),
        auc(roc(dataset$biogrid == "phosphatase", -(dataset$marker))),
        auc(roc(dataset$goldstandard1[!is.na(dataset$goldstandard1)], abs(dataset$marker[!is.na(dataset$goldstandard1)]))),
        auc(roc(dataset$goldstandard2[!is.na(dataset$goldstandard2)], abs(dataset$marker[!is.na(dataset$goldstandard2)]))),
        auc(roc(dataset$goldstandard3[!is.na(dataset$goldstandard3)], abs(dataset$marker[!is.na(dataset$goldstandard3)])))
    )
    paste(c("any", "parca", "fasolo", "biogrid", "kinase", "phosphatase", "gold1", "gold2", "gold3"), aucs, sep="\t", collapse="\n")
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

aucs_P = evaluate_auc_P(P_eval)
aucs_T = evaluate_auc_T(T_eval)
aucs_T_masked = evaluate_auc_T(T_eval_masked)
cat("KP aucs", aucs_P, "TF aucs", aucs_T, "TF maskeed aucs", aucs_T_masked, sep="\n")
cat("KP", evaluate_P(P_eval), sep="\n")
cat("TF", evaluate_T(T_eval), sep="\n")
cat("TF masked", evaluate_T(T_eval_masked), sep="\n")


# plot(roc(P_eval$goldstandard[!is.na(P_eval$goldstandard)], abs(P_eval$marker[!is.na(P_eval$goldstandard)])))




