#!/usr/bin/env Rscript

# packages
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(verification)) # roc auc with p-value

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
# combine pvals using fishers method https://en.wikipedia.org/wiki/Fisher%27s_method 
# chisq = -2 sum(ln(p-values)); pval = 1-pchisq(chisq, df=2length(p-values))
fisher.method.log = function(pvals) {
    df = 2*length(pvals)
    pchisq(-2*sum(log(pvals),na.rm=TRUE),df,lower.tail=FALSE,log.p=TRUE)
}


WP_fnames = commandArgs(trailingOnly=T)
# WP_fnames = "~/cwd/data/inference/02/WP_infer.mat"
P_fname = "~/cwd/data/evaluation/P_eval.tsv"
KP_fname = "~/cwd/data/network/KP.txt"
TF_fname = "~/cwd/data/network/TF.txt"
V_fname = "~/cwd/data/network/V.txt"
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
V = read.vector(V_fname)
PT = c(KP,TF)
nP = length(KP)
nV = length(V)
nT = length(TF)
nO = nV-length(PT)

rundir = getwd()

for (WP_fname in WP_fnames) {
    cat(WP_fname, "\n")
    setwd(dirname(WP_fname))
    WP_fname = basename(WP_fname)
    WP = read.matrix(WP_fname)
    if (all(dim(WP) == nV)) {
        # is square, use the relevant part
        WP = WP[1:(nP+nT),1:nP]
    }
    colnames(WP) = KP
    rownames(WP) = c(KP, TF)
    P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
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
    
    
    # only analysing positives, since we cannot actually rule out negatives
    
    venndata = data.frame(
        potential = TRUE,
        biogrid = P_eval$biogrid != "",
        fasolo = !is.na(P_eval$fasolo),
        parca = !is.na(P_eval$parca),
        fiedler = P_eval$fiedler != "",
        yeastkid = !is.na(P_eval$yeastkid) & P_eval$yeastkid > 4.52,
        ptmod = !is.na(P_eval$ptmod) & P_eval$ptmod > 250,
        ptacek = !is.na(P_eval$ptacek)
    )
    venndata$litterature = venndata$biogrid | venndata$fasolo | venndata$parca | venndata$fiedler
    venndata$curated = venndata$yeastkid | venndata$ptmod
    venndata$known = venndata$litterature | venndata$curated
    venndata$invitro = venndata$known | venndata$ptacek
    
    # index indicating which edge is inferred most strongly to most weakly
    marker_order = order(abs(P_eval$marker), decreasing=T)
    
    p.selection = function(k, selection) {
        intersect = sum(venndata[marker_order[1:k],selection])
        m = sum(venndata[,selection])
        n = nrow(venndata) - m
        list(value=intersect, p=phyper(intersect, m, n, k, lower.tail=FALSE))
    }
    
    quantiles = seq(0.01, 0.333, 0.0025)
    quantiles_plot = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.25, 0.333)
    
    p.selection.q = function(selection, quantiles) {
        out = data.frame()
        for (k in quantiles*nrow(venndata)) {
            out = rbind(out, p.selection(k, selection))
        }
        data.frame(test=paste(selection, quantiles), selection=selection, quantile=quantiles, out)
    }
    
    known.p = p.selection.q("known", quantiles)
    litterature.p = p.selection.q("litterature", quantiles)
    invitro.p = p.selection.q("invitro", quantiles)
    positive.intersections = rbind(known.p, litterature.p, invitro.p)
    best = positive.intersections[which.min(positive.intersections$p),]
    score = -log10(best$p)
    
    eval.table = rbind(cortable, aucs_P, positive.intersections[positive.intersections$quantile%in%quantiles_plot, c("test", "value", "p")])
    eval.table = eval.table[order(eval.table$p),]  # sort by p-value
    
    write.table(eval.table, "evaluation.tsv", sep="\t", quote=F, row.names=F)
    write.table(score, "score.txt", quote=F, row.names=F, col.names=F)
    
    ### Plotting
    selection = as.character(best$selection)
    k = best$quantile * nrow(venndata)
    inferred = unwhich(marker_order[1:k], dim=nrow(P_eval))
    eulerdata = euler(data.frame(venndata[,c("potential", selection)], inferred=inferred), shape="ellipse")
    pdf("euler.pdf")
    print(plot(eulerdata, labels=T, quantities=T))
    dev.off()
    
    fisher.ps = c(fisher.method.log(known.p$p), fisher.method.log(litterature.p$p), fisher.method.log(invitro.p$p))
    pltdf = list(known.p, litterature.p, invitro.p)[[which.min(fisher.ps)]]
    
    second_axis = dup_axis(name="substrates/KP", breaks=quantiles_plot, labels=round(quantiles_plot*nrow(venndata)/length(KP),1))

    plt = ggplot(data=pltdf, aes(x=quantile, y=-log10(p))) +
        geom_line() +
        scale_y_continuous(breaks=seq(0,18,3)) +
        scale_x_continuous(breaks=quantiles_plot, labels=gsub("0\\.","\\.",quantiles_plot), sec.axis=second_axis) +
        theme_linedraw() +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") + 
        annotate("text", x=.31, y=2.01, label="p=0.05", size=3) +
        ylab(expression("-"*log[10]*" p")) +
        theme(panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_line(colour = "lightgray"))
    
    ggsave("hyperp.pdf", plot=plt, width=7, height=3)
    
    
    setwd(rundir)  # go back to so relative dirs for other files still work
}


