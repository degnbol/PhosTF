#!/usr/bin/env Rscript

# packages
library(reshape2)
library(Matrix)
library(eulerr)
# suppressPackageStartupMessages(library(pROC))
library(verification)

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
KP_fname = "~/cwd/data/network/KP.txt"
TF_fname = "~/cwd/data/network/TF.txt"
V_fname = "~/cwd/data/network/V.txt"


P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
WP = read.matrix(WP_fname)
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
V = read.vector(V_fname)
PT = c(KP,TF)
colnames(WP) = KP
rownames(WP) = c(KP, TF)
# swap rownames and colnames columns so the order will be source then target
KP_edges = melt_matrix(WP)[,c(2,1,3)]
colnames(KP_edges) = c("P",  "Target", "marker")
# sanity check
stopifnot(all(as.character(KP_edges$P) == as.character(P_eval$Source)))
stopifnot(all(KP_edges$Target == P_eval$Target))

# hold all data in eval sets
P_eval$marker = KP_edges$marker

# remove diagonals
P_eval = P_eval[as.character(P_eval$Source) != as.character(P_eval$Target),]

evaluate_cor = function(dataset, cor_names, cor_names_pos, cor_names_neg) {
    out = data.frame()
    for (name in c(cor_names, cor_names_pos, cor_names_neg)) {
        valid = !is.na(dataset[name])
        if (name%in%cor_names)   marker = abs(dataset[valid, "marker"])
        if (name%in%cor_names_pos) marker = +(dataset[valid, "marker"])
        if (name%in%cor_names_neg) marker = -(dataset[valid, "marker"])
        PCC = cor.test(dataset[valid,name], marker)[c("estimate", "p.value")]
        out = rbind(out, data.frame(test=paste(name, "cor"), value=PCC[[1]][[1]], p=PCC[[2]]))
    }
    out
}

cor_names = c("yeastkid", "reaction", "ptmod", "expression", "catalysis", "netphorest", "networkin", "networkin STRING", 
              "networkin_biogrid", "undirected", "EMAP", "n_datasets")
cor_names_pos = c("activation")
cor_names_neg = c("inhibition")
cortable = evaluate_cor(P_eval, cor_names, cor_names_pos, cor_names_neg)

aucs = rbind(
    data.frame(roc.area(P_eval$n_datasets > 0, abs(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(!is.na(P_eval$parca), abs(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(!is.na(P_eval$fasolo), abs(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(P_eval$fiedler != "", abs(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(P_eval$biogrid != "", abs(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(P_eval$biogrid == "kinase", +(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(P_eval$biogrid == "phosphatase", -(P_eval$marker))[c("A","p.value")]),
    data.frame(roc.area(P_eval$goldstandard1[!is.na(P_eval$goldstandard1)], abs(P_eval$marker[!is.na(P_eval$goldstandard1)]))[c("A","p.value")]),
    data.frame(roc.area(P_eval$goldstandard2[!is.na(P_eval$goldstandard2)], abs(P_eval$marker[!is.na(P_eval$goldstandard2)]))[c("A","p.value")]),
    data.frame(roc.area(P_eval$goldstandard3[!is.na(P_eval$goldstandard3)], abs(P_eval$marker[!is.na(P_eval$goldstandard3)]))[c("A","p.value")]),
    data.frame(roc.area(P_eval$goldstandard4[!is.na(P_eval$goldstandard4)], abs(P_eval$marker[!is.na(P_eval$goldstandard4)]))[c("A","p.value")])
)
titles = c("any", "parca", "fasolo", "fiedler", "biogrid", "kinase", "phosphatase", "gold1", "gold2", "gold3", "gold4")
aucs_P = cbind(test=paste(titles, "auc"), aucs)
colnames(aucs_P)[2:3] = c("value", "p")

View(rbind(cortable, aucs_P))

# only analysing positives

venndata = data.frame(
    potential = TRUE,
    biogrid = P_eval$biogrid != "",
    fasolo = !is.na(P_eval$fasolo),
    parca = !is.na(P_eval$parca),
    fiedler = P_eval$fiedler != "",
    yeastkid = !is.na(P_eval$yeastkid) & P_eval$yeastkid > 4.52,
    ptmod = !is.na(P_eval$ptmod) & P_eval$ptmod > 250,
)
venndata$litterature = venndata$biogrid | venndata$fasolo | venndata$parca | venndata$fiedler
venndata$curated = venndata$yeastkid | venndata$ptmod
venndata$known = venndata$litterature | venndata$curated


make_euler = function(N, select=c("potential", "known", "inferred")) {
    
    
    euler(venndata[,select], shape="ellipse")
}

inferred = unwhich(order(abs(P_eval$marker), decreasing=T)[1:N], dim=nrow(P_eval))

intersect_3000 = sum(make_euler(3000)$original.values[paste("potential", "known", "inferred", sep="&")])
missed_3000 = sum(make_euler(3000)$original.values[paste("potential", "known", sep="&")])
intersect_5000 = sum(make_euler(5000)$original.values[paste("potential", "known", "inferred", sep="&")])
missed_5000 = sum(make_euler(5000)$original.values[paste("potential", "known", sep="&")])
cat(paste0("euler intersect ", c(3000, 5000), ":\t", c(intersect_3000, intersect_5000), collapse="\n"))


### Plotting

plot(make_euler(5000), labels=T, quantities=T)





