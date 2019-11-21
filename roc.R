#!/usr/bin/env Rscript

# packages
library(ggplot2)
library(plotROC)
library(dplyr)
library(extrafont)
# font_import(); loadfonts() # https://cran.r-project.org/web/packages/extrafont/README.html
# fonts() # see available fonts

# functions
read_vec = function(fname) {as.vector(as.matrix(read.table(fname)))}
get_pos = function(v) {out = +v; out[out <= 0] = 0; out}
get_neg = function(v) {out = -v; out[out <= 0] = 0; out}


main = function(true_fnames, marker_fnames, outfname) {
    # read
    trues = lapply(true_fnames, read_vec)
    markers = lapply(marker_fnames, read_vec)
    Ps  = lapply(trues, abs)
    PKs = lapply(trues, get_pos)
    PPs = lapply(trues, get_neg)
    
    dfs = list()
    for (i in 1:length(trues)) {
        dfs[[i]] = rbind(data.frame(D=Ps[[i]],  M=abs(markers[[i]]), sample=i, type="P"),
                         data.frame(D=PKs[[i]], M=+markers[[i]], sample=i, type="PK"),
                         data.frame(D=PPs[[i]], M=-markers[[i]], sample=i, type="PP"))
    }
    df = bind_rows(dfs)
    df$sample = as.factor(df$sample)
    
    aucs = calc_auc(ggplot(df, aes(d=D, m=M, color=type)) + geom_roc())$AUC
    labels = paste0(c("P  (", "PK (", "PP ("), "AUC=", round(aucs,3), ")")
    
    plot = ggplot(df, aes(d=D, m=M, color=type)) +
        geom_abline(slope=1, intercept=0, color="lightgray") +
        geom_roc(n.cuts=0, aes(fill=sample), linealpha=.2) + 
        geom_roc(n.cuts=0) + 
        style_roc(theme=theme_bw, guide=F, xlab="FPR", ylab="TPR")
    
    # legend and title
    font = "DejaVu Sans Mono"
    plot = plot +
        scale_color_manual(labels=labels, values=c("darkorchid", "red", "blue")) +
        theme(legend.text=element_text(family=font), legend.title=element_blank(), 
              legend.position=c(1,0), legend.justification=c(1,0), legend.box.margin=margin(1, 1, 1, 1)) +
        ggtitle(basename(getwd()))
    
    # plot
    
    ggsave(outfname, plot=plot, width=4.5, height=4.6)
    embed_fonts(outfname)
}


example = function() {
    # filenames
    setwd("/Users/christian/GoogleDrev/PKTFX/data/dream/archive/LB01_100")
    true_fname = "WP.mat"
    marker_fname = "WP_infer.mat"
    true_fnames = list.files(pattern=true_fname, recursive=T)
    marker_fnames = list.files(pattern=marker_fname, recursive=T)
    
    main(true_fnames, marker_fnames)
}

# Example run: roc.R */WP.mat */WP_infer.mat WP.pdf
# first half of infiles are the trues, second half are the markers, except for the last argument that is the output
args = commandArgs(trailingOnly=TRUE)
n.args = length(args)
stopifnot(n.args %% 2 == 1)
true_fnames_idx = 1:((n.args-1)/2)
marker_fnames_idx = ((n.args-1)/2+1):(n.args-1)
outfname_idx = n.args
main(args[true_fnames_idx], args[marker_fnames_idx], args[outfname_idx])
