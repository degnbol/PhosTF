#!/usr/bin/env Rscript

# count number of GO terms in common for each pair of genes

library(data.table)
library(ggplot2)

flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
lower_tri = function(DT) as.matrix(DT)[lower.tri(DT)]
melt_ = function(x) melt(x, id.vars="ORF", variable.name="KP", value.name="w")
wilcox_p = function(dataset) {
    wilcox.test(dataset[infer==TRUE,shared], dataset$shared, alternative="g")$p.value
}

# constants
KP_color = "#bd61b6"
TF_color = "#75b42f"

KPs = read.vector("~/cwd/data/network/KP.txt")
TFs = read.vector("~/cwd/data/network/TF.txt")
Vs = read.vector("~/cwd/data/network/V.txt")
PTs = c(KPs, TFs)
shared_GO = fread("shared_GO.tsv", sep="\t")


# WP_fnames = "~/cwd/data/inference/01/WP_infer.mat"  # test
WP_fnames = commandArgs(trailingOnly=TRUE)
n_top_KP = length(KPs) * 100
n_top_TF = length(TFs) * 100

init_dir = getwd()
for(WP_fname in WP_fnames) {
    setwd(dirname(WP_fname))
    WP_fname = basename(WP_fname)
    
    WP = fread(WP_fname, sep=" ", col.names=KPs)
    WP = data.table(ORF=PTs, abs(WP), key="ORF")
    infers = rbind(melt_(WP[KPs,])[order(-w)[1:n_top_KP],], melt_(WP[TFs,])[order(-w)[1:n_top_TF],])
    infers[,w:=NULL][,infer:=T]
    
    shared_GO[infers,on=c("ORF","KP"),infer:=infer]
    shared_GO$infer[is.na(shared_GO$infer)] = F
    
    p.values = data.frame(Substrate=c("KP","TF"),
        p=c(wilcox_p(shared_GO[Aspect=="P" & Substrate == "KP",]),
            wilcox_p(shared_GO[Aspect=="P" & Substrate == "TF",])))
    fwrite(list(p.values$p), "shared_GO_p.tsv", sep="\t")
    
    # mean shared GO terms for different interesting groupings
    shared_GO_dens = shared_GO[,.N,by=c("infer", "Aspect", "Substrate", "shared")]
    shared_GO_dens[shared_GO[,.(out_of=.N),by=c("infer", "Aspect", "Substrate")],
                   on=c("infer", "Aspect", "Substrate"), out_of:=out_of]
    shared_GO_dens[,density:=N/out_of]
    
    plt = ggplot(shared_GO_dens[Aspect=="P",]) +
        geom_col(mapping=aes(x=shared, y=density, fill=infer), position="dodge", width=.8) +
        theme_linedraw() +
        scale_fill_manual(values=c("darkgray", "black"), breaks=list(TRUE,FALSE), labels=c("inferred", "potential")) +
        scale_x_continuous(breaks=c(0,1,2,5,10), limits=c(-.5,14.5), expand=c(0,0)) +
        facet_grid(vars(Substrate)) + 
        xlab("shared GO pathways") + 
        theme(legend.title=element_blank(), 
              panel.grid.major=element_line(colour="gray"), 
              panel.grid.minor=element_line(colour="lightgray")) +
        geom_text(data=p.values, mapping=aes(label=sprintf("p=%.3g",p)), x=12.5, y=.7)
    
    ggsave("shared_GO.pdf", plot=plt, width=6, height=3)
    
    setwd(init_dir)
}











