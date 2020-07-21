#!/usr/bin/env Rscript

library(data.table)
library(fdrtool)
library(reshape2)
library(ggplot2)

## functions

flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
read.matrix = function(x) as.matrix(read.table(x))
melt_matrix = function(x) {
    out = melt(as.matrix(x))
    colnames(out) = c("rownames", "colnames", "value")
    out
}



# WP_fname = "~/cwd/data/inference/74/WP_infer_1.mat"
WP_fnames = commandArgs(trailingOnly=T)
KP = fread("~/cwd/data/network/KP_protein.tsv")$ORF
TF = read.vector("~/cwd/data/network/TF.txt")
V = read.vector("~/cwd/data/network/V_protein.txt")
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
    # swap rownames and colnames columns so the order will be source then target
    KP_edges = melt_matrix(WP)[,c(2,1,3)]
    colnames(KP_edges) = c("KP", "Target", "marker")
    
    # remove diagonals
    KP_edges = KP_edges[as.character(KP_edges$KP) != as.character(KP_edges$Target),]
    
    
    KP2KP.idx = KP_edges$Target%in%KP
    KP2TF.idx = KP_edges$Target%in%TF
    stopifnot(all(KP2KP.idx | KP2TF.idx))
    stopifnot(!any(KP2KP.idx & KP2TF.idx))
    
    
    KP_edges$q = NA
    KP_edges$q[KP2KP.idx] = fdrtool(KP_edges$marker[KP2KP.idx], plot=FALSE)$qval
    KP_edges$q[KP2TF.idx] = fdrtool(KP_edges$marker[KP2TF.idx], plot=FALSE)$qval
    
    write.table(KP_edges, "KP_edges.tsv", sep="\t", quote=F, row.names=F)
    KP_edges$infer = KP_edges$q < .05
    cat(sum(KP_edges$infer[KP2KP.idx]), " ", sum(KP_edges$infer[KP2TF.idx]), "\n")
    
    plt = ggplot(KP_edges, aes(marker, fill=infer)) + 
        geom_histogram(binwidth=.001) +
        theme_linedraw() +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        scale_color_manual(values=c("transparent", "black")) +
        theme(panel.grid.major.x=element_line(colour="gray"), panel.grid.minor.x=element_line(colour="lightgray"))
    
    ggsave("WPFDR.pdf", plot=plt)
    
    setwd(rundir)
}



