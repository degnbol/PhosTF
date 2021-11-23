#!/usr/bin/env Rscript
# packages
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(extrafont))
# font_import(); loadfonts() # https://cran.r-project.org/web/packages/extrafont/README.html
# fonts() # see available fonts
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(here))
library(Cairo)

# functions
read = function(fname) {as.matrix(read.table(fname))}
read_vec = function(fname) {as.vector(as.matrix(read.table(fname)))}
get_pos = function(v) {out = +v; out[out <= 0] = 0; out}
get_neg = function(v) {out = -v; out[out <= 0] = 0; out}


main = function(true_fnames, marker_fnames, tf_true_fnames, tf_marker_fnames, outfname, title) {
    # read
    trues = lapply(true_fnames, read)
    tf_trues = lapply(tf_true_fnames, read)
    markers = lapply(marker_fnames, read)
    tf_markers = lapply(tf_marker_fnames, read)
    Ps  = lapply(trues, abs)
    Ts  = lapply(tf_trues, abs)
    PKs = lapply(trues, get_pos)
    PPs = lapply(trues, get_neg)
    
    dfs = list()
    for (i in 1:length(trues)) {
        nP  = ncol(trues[[i]])
        nPT = nrow(trues[[i]])
        n   = nrow(tf_trues[[i]])
        M_abs = abs(markers[[i]])
        M_tf  = abs(tf_markers[[i]])
        M_pos = +markers[[i]]
        M_neg = -markers[[i]]
        P  = Ps[[i]]
        T_ = Ts[[i]]
        PK = PKs[[i]]
        PP = PPs[[i]]
        dfs[[i]] = rbind(
            data.frame(D=as.vector(P[1:nP,]),        M=as.vector(M_abs[1:nP,]),       sample=i, source="KP ",  target="KP"),
            data.frame(D=as.vector(PK[1:nP,]),       M=as.vector(M_pos[1:nP,]),       sample=i, source="KP+",  target="KP"),
            data.frame(D=as.vector(PP[1:nP,]),       M=as.vector(M_neg[1:nP,]),       sample=i, source="KP-",  target="KP"),
            data.frame(D=as.vector(P[(nP+1):nPT,]),  M=as.vector(M_abs[(nP+1):nPT,]), sample=i, source="KP ",  target="TF"),
            data.frame(D=as.vector(PK[(nP+1):nPT,]), M=as.vector(M_pos[(nP+1):nPT,]), sample=i, source="KP+",  target="TF"),
            data.frame(D=as.vector(PP[(nP+1):nPT,]), M=as.vector(M_neg[(nP+1):nPT,]), sample=i, source="KP-",  target="TF"),
            data.frame(D=as.vector(T_),              M=as.vector(M_tf),               sample=i, source="TF ",  target="V "))
    }
    df = bind_rows(dfs)
    df$sample = as.factor(df$sample)
    lvls = c("TF →V ", "KP →KP", "KP+→KP", "KP-→KP", "KP →TF", "KP+→TF", "KP-→TF")
    cols = c("#005000", "#B663B1", "#B663B1", "#B663B1", "#6DB845", "#6DB845", "#6DB845")
    # "31" means 3 units (linewidth) solid, then 1 unit gap. This is the way to control line width separately from size, i.e. how far between dots etc.
    ltyp = c("solid", "solid", "31", "11", "solid", "31", "11")
    df$edge = factor(paste0(df$source, "→", df$target), levels=lvls)
    
    aucs = calc_auc(ggplot(df, aes(d=D, m=M, color=edge)) + geom_roc())$AUC
    labels = paste0(levels(df$edge), "   ", sprintf("%.3f", aucs), " ")
    
    plot = ggplot(df, aes(d=D, m=M, color=edge, linetype=edge)) +
        geom_abline(slope=1, intercept=0, color="lightgray") +
        geom_roc(n.cuts=0, size=.8) +
        style_roc(theme=theme_bw, guide=F, xlab="False Positive Rate", ylab="True Positive Rate")
    # + geom_roc(n.cuts=0, aes(fill=sample), linealpha=.1) # transparent individual curves
    
    # legend and title
    font = "DejaVu Sans Mono"
    legend_title = "   Src.→Trg.  AUC"
    plot = plot +
        scale_color_manual(labels=labels, values=cols, name=legend_title) +
        scale_linetype_manual(labels=labels, values=ltyp, name=legend_title) +
        theme(legend.text=element_text(family=font, size=9), legend.title=element_text(family=font, size=9), legend.position=c(1,0), 
              legend.justification=c(1,0), legend.box.margin=margin(0, 1, 1, 0), legend.margin=margin(0,0,0,0)) +
        ggtitle(title)
    
    # plot
    
    ggsave(outfname, plot=plot, width=4.2, height=4.28, device=cairo_pdf)
    embed_fonts(outfname)
}


example = function() {
    # filenames
    setwd(paste0(here(), "/data/dream/archive/LB01_100"))
    true_fname = "WP.mat"
    tf_true_fname = "WT.mat"
    marker_fname = "WP_infer.mat"
    tf_marker_fname = "WT_infer.mat"
    true_fnames = list.files(pattern=true_fname, recursive=T)
    tf_true_fnames = list.files(pattern=tf_true_fname, recursive=T)
    marker_fnames = list.files(pattern=marker_fname, recursive=T)
    tf_marker_fnames = list.files(pattern=tf_marker_fname, recursive=T)
    
    main(true_fnames, marker_fnames, tf_true_fnames, tf_marker_fnames, "ROC.pdf", "ROC")
}

get_parser = function() {
    parser = ArgumentParser(description="Plot a ROC curve.")
    parser$add_argument("-t", "--true",   type="character", nargs="+", help="Matrix files with true values. Same number of files as --marker.")
    parser$add_argument("-m", "--marker", type="character", nargs="+", help="Matrix files with marker values. Same number of files as --true.")
    parser$add_argument("--tf-true",   type="character", nargs="*", help="Matrix files with true values for TF. Same number of files as --tf-marker.")
    parser$add_argument("--tf-marker", type="character", nargs="*", help="Matrix files with marker values for TF. Same number of files as --tf-true.")
    parser$add_argument("-o", "--out", type="character", default="ROC.pdf", help="output file name for roc plot. (default = %(default)s)")
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
    if(length(args$tf_true) != length(args$tf_marker)) {
        stop("The same number of --tf-true and --tf-marker files must be supplied.")
    }
    args
}

args = get_args()

main(args$true, args$marker, args$tf_true, args$tf_marker, args$out, args$title)

