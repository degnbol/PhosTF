#!/usr/bin/env Rscript

library(fdrtool)


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
KP = read.vector("~/cwd/data/network/KP.txt")
TF = read.vector("~/cwd/data/network/TF.txt")
V = read.vector("~/cwd/data/network/V.txt")
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
 
    
    
    KP_edges$q = fdrtool(abs(KP_edges$marker), plot=FALSE)$qval
    cat("raw < .1", sum(fdrtool(KP_edges$marker, plot=FALSE)$qval < .05) / 199, "\n")
    cat("abs < .2", sum(KP_edges$q < .2) / 199, "\n")
    KP_edges$infer = KP_edges$q < .2
    
    
    plt = ggplot(KP_edges, aes(marker, fill=infer)) + 
        geom_histogram(binwidth=.005) +
        theme_linedraw() +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        scale_color_manual(values=c("transparent", "black")) +
        theme(panel.grid.major.x=element_line(colour="gray"), panel.grid.minor.x=element_line(colour="lightgray"))
    
    ggsave("WPFDR.pdf", plot=plt)
    
    setwd(rundir)
}



