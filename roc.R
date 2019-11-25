#!/usr/bin/env Rscript
# packages
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(extrafont))
# font_import(); loadfonts() # https://cran.r-project.org/web/packages/extrafont/README.html
# fonts() # see available fonts
suppressPackageStartupMessages(library(argparse))

# functions
read_vec = function(fname) {as.vector(as.matrix(read.table(fname)))}
get_pos = function(v) {out = +v; out[out <= 0] = 0; out}
get_neg = function(v) {out = -v; out[out <= 0] = 0; out}


main = function(true_fnames, marker_fnames, outfname, title) {
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
        geom_roc(n.cuts=0) +
        style_roc(theme=theme_bw, guide=F, xlab="FPR", ylab="TPR")
    # + geom_roc(n.cuts=0, aes(fill=sample), linealpha=.1) # transparent individual curves
    
    # legend and title
    font = "DejaVu Sans Mono"
    plot = plot +
        scale_color_manual(labels=labels, values=c("darkorchid", "red", "blue")) +
        theme(legend.text=element_text(family=font), legend.title=element_blank(), 
              legend.position=c(1,0), legend.justification=c(1,0), legend.box.margin=margin(1, 1, 1, 1)) +
        ggtitle(title)
    
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
    
    main(true_fnames, marker_fnames, "ROC.pdf", "ROC")
}


get_parser = function() {
    parser = ArgumentParser(description="Plot a ROC curve.")
    parser$add_argument("-t", "--true",   type="character", nargs="+", help="Matrix files with true values. Same number of files as --marker.")
    parser$add_argument("-m", "--marker", type="character", nargs="+", help="Matrix files with marker values. Same number of files as --true.")
    parser$add_argument("-o", "--out", type="character", default="roc.pdf", help="output file name for roc plot. (default = %(default)s)")
    parser$add_argument("--title", type="character", default=basename(getwd()), help="Title of plot. (default = %(default)s)")
    parser
}

get_args = function() {
    args = get_parser()$parse_args()
    if(is.null(args$true) || is.null(args$marker)) {
        stop("At least one --true and --marker file must be supplied.")
    }
    if(length(args$true) != length(args$marker)) {
        stop("The same number of --true and --marker files must be supplied.")
    }
    args
}

args = get_args()

main(args$true, args$marker, args$out, args$title)
