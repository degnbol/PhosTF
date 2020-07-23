#!/usr/bin/env Rscript

# packages
library(reshape2)
library(Matrix)
library(eulerr)
library(ggplot2)
library(ggpubr)
library(data.table)

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






## settings
# assuming we are masking KP edges based on known phos site from BioGrid
masking_KP = TRUE


# read evaluation files
KP = fread("~/cwd/data/network/KP_protein.tsv")$ORF
TF = read.vector("~/cwd/data/network/TF.txt")
V  = read.vector("~/cwd/data/network/V_protein.txt")
KPTF = c(KP,TF)
nKP = length(KP)
nV = length(V)
nTF = length(TF)
nO = nV-length(KPTF)
KP_eval = fread("~/cwd/data/evaluation/P_eval.tsv", 
                select=c("Source", "Target", "biogrid", "fasolo", "parca", "fiedler", "yeastkid", "ptmod", "ptacek"))
# remove diagonals
KP_eval = KP_eval[Source!=Target,]
# make evaluations Boolean
KP_eval[,biogrid:=biogrid!=""]
KP_eval[,fiedler:=fiedler!=""]
KP_eval[,fasolo:=!is.na(fasolo)]
KP_eval[,parca:=!is.na(parca)]
KP_eval[,ptacek:=!is.na(ptacek)]
# 4.52 corresponds to p<0.05 according to yeastKID
KP_eval[,yeastkid:=!is.na(yeastkid) & yeastkid>4.52]
# string DB scores are confidence from 0 to 1 written as e.g. 800 to mean 0.8.
# higher confidence score is better. 0.5 means that there is 50% prob for false positive. https://string-db.org/cgi/info.pl
# if we use False Discovery Rate = .2 it only gives us a meager 4 edges, all of which are already in other datasets.
# # so we use a dangerous FDR = .4, which gives us 7 edges not already found in another evaluation set.
KP_eval[,ptmod:=!is.na(ptmod) & ptmod>600]
KP_eval[,literature:=biogrid|fasolo|parca|fiedler]
KP_eval[,curated:=yeastkid|ptmod]
KP_eval[,known:=literature|curated]
KP_eval[,invitro:=known|ptacek]
# not evaluating on uninferrable edges
KP_eval$possible = T
if(masking_KP) {
    KP_targets_noknownsite = read.vector("~/cwd/data/network/KP_targets_noknownsite.txt")
    KP_eval[Target%in%KP_targets_noknownsite,possible:=F]
}


#KP_edge_fname = "~/cwd/data/inference/03/KP_edges.tsv"
KP_edge_fnames = commandArgs(trailingOnly=T)
rundir = getwd()

for(evalname in c("known", "literature", "invitro")) {
    plts = list()
    for (i_file in 1:length(KP_edge_fnames)) {
        KP_edge_fname = KP_edge_fnames[i_file]
        cat(KP_edge_fname, "\n")
        setwd(dirname(KP_edge_fname))
        KP_edges = fread(basename(KP_edge_fname), drop="marker")
        KP_edges$infer = KP_edges$q < .05
        KP_edges[,q:=NULL]
        
        # all targets should be protein KP or TF
        stopifnot(all(KP_edges$KP%in%KP))
        stopifnot(all(KP_edges$Target%in%KPTF))
        stopifnot(all(KP_eval$Source%in%KP))
        stopifnot(all(KP_eval$Target%in%KPTF))
        # sanity check
        stopifnot(all(KP_edges$KP == KP_eval$Source))
        stopifnot(all(KP_edges$Target == KP_eval$Target))
        # add eval data
        KP_edges = cbind(KP_edges, KP_eval[,!c("Source", "Target")])
        # filter for possible edges and remove then redundant column
        KP_edges = KP_edges[possible==T,][,possible:=NULL]
        
        KP_edges[Target%in%KP,Target:="KP"]
        KP_edges[Target%in%TF,Target:="TF"]
        KP_edges[,KP:=NULL]
            
        DT = KP_edges[, c("Target", "infer", evalname), with=F]
        ps = c(get_p(DT[Target=="KP"][[evalname]], DT[Target=="KP"][["infer"]]), 
               get_p(DT[Target=="TF"][[evalname]], DT[Target=="TF"][["infer"]]))
        
        setnames(DT, "Target", "substrate")
        setnames(DT, evalname, "eval")
        plts[[i_file]] = plot_square_euler(DT, ps)
            
        setwd(rundir)  # go back to so relative dirs for other files still work
    }
    plt = ggarrange(plotlist=plts, ncol=1, labels=KP_edge_fnames)
    ggsave(paste0("square_euler_", evalname, ".pdf"), plot=plt, width=6.5, height=1.35)
}
