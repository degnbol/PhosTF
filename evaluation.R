#!/usr/bin/env Rscript

# packages
library(reshape2)
library(Matrix)
library(eulerr)
library(lmPerm)
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
unwhich = function(which, dim=max(which)) {
    y = array(logical(length(which)), dim=dim)
    y[which] = TRUE
    y
}

example = TRUE
if (example) {
    setwd("~/cwd/data/inference/02")
    WP_fname = "WP_infer.mat"
} else {
    args = commandArgs(trailingOnly=T)
    WP_fname = args[1]
}

P_fname = "~/cwd/data/evaluation/P_eval.tsv"
T_fname = "~/cwd/data/evaluation/T_eval.tsv"
WT_mask_fname = "~/cwd/data/network/WT_mask.csv"
KP_fname = "~/cwd/data/network/KP.txt"
TF_fname = "~/cwd/data/network/TF.txt"
V_fname = "~/cwd/data/network/V.txt"


P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
T_eval = read.table(T_fname, header=T, sep="\t", quote="", check.names=F)
WP = read.matrix(WP_fname)
WT = read.matrix(WT_fname)
WT_mask = as.matrix(read.table(WT_mask_fname, sep=",", check.names=F))
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
V = read.vector(V_fname)
PT = c(KP,TF)
colnames(WP) = KP
colnames(WT) = TF
stopifnot(all(colnames(WT_mask) == TF))
rownames(WP) = c(KP, TF)
rownames(WT) = V
stopifnot(all(rownames(WT_mask) == V))
# swap rownames and colnames columns so the order will be source then target
KP_edges = melt_matrix(WP)[,c(2,1,3)]
TF_edges = melt_matrix(WT)[,c(2,1,3)]
TF_mask_edges = melt_matrix(WT_mask)[,c(2,1,3)]
colnames(KP_edges) = c("P",  "Target", "marker")
colnames(TF_edges) = c("TF", "Target", "marker")
colnames(TF_mask_edges) = c("TF", "Target", "mask")
# sanity check
all(as.character(KP_edges$P) == as.character(P_eval$Source))
all(KP_edges$Target == P_eval$Target)

# hold all data in eval sets
P_eval$marker = KP_edges$marker

# remove diagonals
P_eval = P_eval[as.character(P_eval$Source) != as.character(P_eval$Target),]


evaluate_cor = function(dataset, cor_names, cor_names_pos, cor_names_neg) {
    out = c()
    for (name in cor_names) {
        valid = !is.na(dataset[name])
        PCC = cor.test(dataset[valid,name], abs(dataset[valid, "marker"]))[c("estimate", "p.value")]
        PCC = summary(lm(dataset[valid,name] ~ abs(dataset[valid, "marker"])))$coefficients[2,3]
        out = c(out, sprintf(paste0(name, "%.4f", sep="\t"), PCC))
    }
    for (name in cor_names_pos) {
        valid = !is.na(dataset[name])
        PCC = lm(dataset[valid,name], +(dataset[valid, "marker"]))$coefficients[2,3]
        out = c(out, sprintf(paste0(name, "%.4f", sep="\t"), PCC))
    }
    for (name in cor_names_neg) {
        valid = !is.na(dataset[name])
        PCC = lm(dataset[valid,name], -(dataset[valid, "marker"]))$coefficients[2,3]
        out = c(out, sprintf(paste0(name, "%.4f", sep="\t"), PCC))
    }
    paste(out, collapse="\n")
}


evaluate_P = function(dataset) {
    cor_names = c("yeastkid", "reaction", "ptmod", "expression", "catalysis", "netphorest", "networkin", "networkin STRING", 
                  "networkin_biogrid", "undirected", "EMAP", "n_datasets")
    cor_names_pos = c("activation")
    cor_names_neg = c("inhibition")
    evaluate_cor(dataset, cor_names, cor_names_pos, cor_names_neg)
}


evaluate_auc_P = function(dataset) {
    aucs = c(
        auc(roc(dataset$n_datasets > 0, abs(dataset$marker), direction="<")),
        auc(roc(!is.na(dataset$parca), abs(dataset$marker), direction="<")),
        auc(roc(!is.na(dataset$fasolo), abs(dataset$marker), direction="<")),
        auc(roc(dataset$fiedler != "", abs(dataset$marker), direction="<")),
        auc(roc(dataset$biogrid != "", abs(dataset$marker), direction="<")),
        auc(roc(dataset$biogrid == "kinase", +(dataset$marker), direction="<")),
        auc(roc(dataset$biogrid == "phosphatase", -(dataset$marker), direction="<")),
        auc(roc(dataset$goldstandard1[!is.na(dataset$goldstandard1)], abs(dataset$marker[!is.na(dataset$goldstandard1)]), direction="<")),
        auc(roc(dataset$goldstandard2[!is.na(dataset$goldstandard2)], abs(dataset$marker[!is.na(dataset$goldstandard2)]), direction="<")),
        auc(roc(dataset$goldstandard3[!is.na(dataset$goldstandard3)], abs(dataset$marker[!is.na(dataset$goldstandard3)]), direction="<")),
        auc(roc(dataset$goldstandard4[!is.na(dataset$goldstandard4)], abs(dataset$marker[!is.na(dataset$goldstandard4)]), direction="<"))
    )
    paste(c("any", "parca", "fasolo", "fiedler", "biogrid", "kinase", "phosphatase", "gold1", "gold2", "gold3", "gold4"), aucs, sep="\t", collapse="\n")
}

aucs_P = evaluate_auc_P(P_eval)
cat("KP aucs", aucs_P, sep="\n")
cat("KP cor", evaluate_P(P_eval), sep="\n")

make_euler = function(N, select=c("potential", "known", "inferred")) {
    venndata = data.frame(
        potential = T,
        biogrid = P_eval$biogrid != "",
        fasolo = !is.na(P_eval$fasolo),
        parca = !is.na(P_eval$parca),
        fiedler = P_eval$fiedler != "",
        yeastkid = !is.na(P_eval$yeastkid) & P_eval$yeastkid > 4.52,
        ptmod = !is.na(P_eval$ptmod) & P_eval$ptmod > 250,
        # netphorest = unwhich(order(P_eval$netphorest, decreasing=T)[1:N], dim=nrow(P_eval))
        # networkin = unwhich(order(P_eval$networkin, decreasing=T)[1:N], dim=nrow(P_eval)),
        # random = sample(c(rep(T,N),rep(F,nrow(P_eval)-N))),
        inferred = unwhich(order(abs(P_eval$marker), decreasing=T)[1:N], dim=nrow(P_eval))
    )
    venndata$litterature = venndata$biogrid | venndata$fasolo | venndata$parca | venndata$fiedler
    venndata$curated = venndata$yeastkid | venndata$ptmod
    venndata$known = venndata$litterature | venndata$curated
    
    euler(venndata[,select], shape="ellipse")
}

intersect_3000 = sum(make_euler(3000)$original.values[paste("potential", "known", "inferred", sep="&")])
missed_3000 = sum(make_euler(3000)$original.values[paste("potential", "known", sep="&")])
intersect_5000 = sum(make_euler(5000)$original.values[paste("potential", "known", "inferred", sep="&")])
missed_5000 = sum(make_euler(5000)$original.values[paste("potential", "known", sep="&")])
cat(paste0("euler intersect ", c(3000, 5000), ":\t", c(intersect_3000, intersect_5000), collapse="\n"))

plot(make_euler(5000), labels=T, quantities=T)





