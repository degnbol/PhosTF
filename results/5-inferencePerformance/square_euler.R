#!/usr/bin/env Rscript
library(Matrix)
library(eulerr)
library(ggplot2)
library(ggpubr)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))
cwd = function(s) paste0(here(), "/", s)

# args
# USE: ./square_euler.R WP_infer*.tsv
# WP_fnames = cwd("results/4-yeastNetworkReconstruction/infer/WP_infer-inner-strict-0.0-0.1.tsv")
WP_fnames = commandArgs(trailingOnly=T)

# functions
flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(cwd(x)))
read.matrix = function(x) as.matrix(read.table(cwd(x)))
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
# evalset is bool vector indicating true and false edges.
# drawn is a bool vector of same size indicating the edges that are inferred.
get_p = function(evalset, drawn) {
    q = sum(evalset[drawn])  # number of correct inferences
    m = sum(evalset) # number of true edges
    n = length(evalset) - m  # number of black balls to pick/"false" edges/potential edges not in set of known
    k = sum(drawn)  # number of inferred edges
    phyper(q, m, n, k, lower.tail=F)
}
# DT has columns substrate, infer, eval. Substrate is "KP" or "TF". infer and eval are Boolean.
# ps is p-value for KP2KP and KP2TF
plot_square_euler = function(DT, ps) {
    ## constants
    KP_color = "#bd61b6"
    TF_color = "#75b42f"
    bw = .8
    hn = .05  # horizontal nudge
    fgscl = 4  # how wide foreground is compared to background
    # version where we nudge foregrounds down 1/4
    fgdg = bw/4
    
    
    ## make background and foreground bars
    bg = DT[       , .N, by=c("substrate", "infer")]
    fg = DT[eval==T, .N, by=c("substrate", "infer")]
    # labels
    bg[,label:=fifelse(infer, "Inferred", "Not inferred")]
    fg[,label:=substrate]
    # bar heights
    bg[, y:=fifelse(infer, -N, N)]
    fg[, y:=fifelse(infer, -N, N)]
    # horizontal adjustment
    bg[, hjust:=fifelse(infer, -hn, 1+hn)]
    fg[, hjust:=fifelse(infer, 1+hn, -hn)]
    # x
    bg[, x:=fifelse(substrate == "KP", 1, 2) + bw/2-.01]
    fg[, x:=fifelse(substrate == "KP", 1, 2) - fgdg]
    # colors
    fg[, color:=fifelse(substrate == "KP", KP_color, TF_color)]
    
    ggplot(mapping=aes(x=substrate, y=y)) + 
        # outline
        geom_bar(data=bg, mapping=aes(fill=label), stat="identity", position="identity", width=bw, fill="transparent", color="black", size=2) +
        # gray and white background bars
        geom_bar(data=bg, mapping=aes(fill=label), stat="identity", position="identity", width=bw) +
        # colored foreground bars with eval data
        geom_bar(data=fg, mapping=aes(fill=substrate, y=y*fgscl), stat="identity", position=position_nudge(x=-fgdg), width=bw/fgscl, alpha=.5) +
        scale_fill_manual(values=list(KP=KP_color, TF=TF_color, Inferred="lightgray", `Not inferred`="white")) +
        # labels saying "Inferred" and "Not inferred"
        annotate("text", x=bg$x, vjust=1,   y=bg$y, hjust=bg$hjust, label=bg$label, fontface="italic", alpha=.7) +
        # numbers under those labels
        annotate("text", x=bg$x, vjust=2.2, y=bg$y, hjust=bg$hjust, label=bg$N,     fontface="bold", alpha=.8) +
        # labels "Known"
        annotate("text", x=c(1,2)-fgdg+bw/fgscl, vjust=0, hjust=.5, y=0, label="Known", fontface="italic", color=c(KP_color, TF_color)) +
        # labels with counts
        annotate("text", x=fg$x, vjust=0.5, y=fg$y*fgscl, hjust=fg$hjust, label=fg$N, fontface="bold", color=fg$color) +
        # p-values
        annotate("text", x=c(1,2)-fgdg, vjust=0.5, hjust=0, y=6000, label=sprintf("p=%.3g",ps), color=c(KP_color, TF_color)) +
        coord_flip() + theme_void() + 
        theme(legend.position="none", axis.text.y=element_text(), plot.margin=margin(t=0,b=0,l=0,r=0,unit="pt")) +
        scale_x_discrete(labels=c("d(KP,KP)", "d(KP,TF)"))
}

# read evaluation files
KP_eval = fread(cwd("results/4-yeastNetworkReconstruction/P_edges/P_edges.tsv"))

plts = list()
for (i_file in 1:length(WP_fnames)) {
    cat(basename(WP_fnames[i_file]), "\t")
    WP = fread(WP_fnames[i_file])
    KPs = names(WP)[2:ncol(WP)]
    # melt
    names(WP)[1] = "Target"
    KP_edges = melt(WP, variable.name="P", id.vars="Target", value.name="pred")
    KP_edges = KP_edges[P!=Target]
    KP_edges[, substrate:="TF"]
    KP_edges[Target%in%KPs, substrate:="KP"]
    
    # join with eval data
    KP_eval[, eval:=TRUE]
    DT = KP_eval[KP_edges, on=c("P", "Target")]
    DT[is.na(eval), eval:=FALSE]
    
    ps = c(1, 1)
    for (thres in c(.95, .9, .8)) {
        DT[, infer:=abs(pred) > quantile(abs(pred), thres), by=substrate]
        ps = pmin(ps, c(get_p(DT[substrate=="KP", eval], DT[substrate=="KP", infer]), 
                        get_p(DT[substrate=="TF", eval], DT[substrate=="TF", infer])))
        DT[, infer:=abs(pred) > quantile(abs(pred), thres)]
        ps = pmin(ps, c(get_p(DT[substrate=="KP", eval], DT[substrate=="KP", infer]), 
                        get_p(DT[substrate=="TF", eval], DT[substrate=="TF", infer])))
    }
        
    
    # plts[[i_file]] = plot_square_euler(DT, ps)
    cat(paste(sprintf("%.3f", -log10(ps)), collapse=" "), "\n")
    
}

# plt = ggarrange(plotlist=plts, ncol=1, labels=WP_fnames)
# outfname = paste0(dirname(WP_fnames[i_file]), "/square_euler.pdf")
# ggsave(outfname, plot=plt, width=6.5, height=1.35*length(plts))

