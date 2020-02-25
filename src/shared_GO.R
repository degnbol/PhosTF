#!/usr/bin/env Rscript

# count number of GO terms in common for each pair of genes

library(data.table)
library(ggplot2)

flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
lower_tri = function(DT) as.matrix(DT)[lower.tri(DT)]
melt_ = function(x) melt(x, id.vars="ORF", variable.name="KP", value.name="w")
wilcox_p = function(dataset) {
    wilcox.test(dataset[infer==TRUE,shared], dataset[infer==FALSE,shared], alternative="g")$p.value
}

# constants
KP_color = "#bd61b6"
TF_color = "#75b42f"

KPs = read.vector("~/cwd/data/network/KP.txt")
TFs = read.vector("~/cwd/data/network/TF.txt")
Vs = read.vector("~/cwd/data/network/V.txt")
PTs = c(KPs, TFs)
nKP = length(KPs)
nTF = length(TFs)
nPT = nKP+nTF


KP_edges_fname = "~/cwd/data/inference/74/KP_edges.tsv"
# KP_edges_fnames = commandArgs(trailingOnly=TRUE)

init_dir = getwd()
for(KP_edges_fname in KP_edges_fnames) {
    cat(KP_edges_fname, "\n")
    setwd(dirname(KP_edges_fname))
    KP_edges = fread(basename(KP_edges_fname), sep="\t")
    KP_edges[,infer:=q<.05]
    colnames(KP_edges)[colnames(KP_edges)=="Target"] = "ORF"
    shared_GO = fread("~/cwd/data/go/shared_GO.tsv", sep="\t")
    
    shared_GO[KP_edges,on=c("ORF","KP"),infer:=infer]
    stopifnot(!any(is.na(shared_GO$infer)))
    
    p.values = data.frame(Substrate=c("KP","TF"),
                          p=c(wilcox_p(shared_GO[Aspect=="P" & Substrate == "KP",]),
                              wilcox_p(shared_GO[Aspect=="P" & Substrate == "TF",])))
    
    # mean shared GO terms for different interesting groupings
    shared_GO_dens = shared_GO[,.N,by=c("infer", "Aspect", "Substrate", "shared")]
    shared_GO_dens[shared_GO[,.(out_of=.N),by=c("infer", "Aspect", "Substrate")],
                   on=c("infer", "Aspect", "Substrate"), out_of:=out_of]
    shared_GO_dens[,density:=N/out_of]
    
    plt = ggplot(shared_GO_dens[Aspect=="P",]) +
        geom_col(mapping=aes(x=shared, y=density, fill=infer), position="dodge", width=.8) +
        theme_linedraw() +
        scale_fill_manual(name="inferred", values=c("darkgray", "black"), breaks=list(TRUE,FALSE), labels=c("yes", "no")) +
        scale_x_continuous(breaks=c(0,1,2,5,10), limits=c(-.5,14.5), expand=c(0,0), minor_breaks=c(3,4,6,7,8,9)) +
        scale_y_log10() +
        facet_grid(vars(Substrate)) + 
        xlab("shared GO pathways") + 
        theme(panel.grid.major=element_line(colour="gray"), 
              panel.grid.minor=element_line(colour="lightgray")) +
        geom_text(data=p.values, mapping=aes(label=sprintf("p=%.3g",p)), x=11.5, y=.7)
    
    ggsave("shared_GO.pdf", plot=plt, width=6, height=3)
    
    setwd(init_dir)
}











