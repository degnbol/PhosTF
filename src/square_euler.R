#!/usr/bin/env Rscript

# packages
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(ggplot2))
library(latex2exp)

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

# convert tex code to expression object (use \\ for \ and combine expressions with list() instead of c())
tex = function(x) unname(TeX(paste0("$",x)))


# KP_edge_fname = "~/cwd/data/inference/74/KP_edges.tsv"
KP_edge_fnames = commandArgs(trailingOnly=T)
P_fname = "~/cwd/data/evaluation/P_eval.tsv"
P_eval_noknownsite_fname = "~/cwd/data/evaluation/KP_targets_noknownsite.txt"
KP_fname = "~/cwd/data/network/KP.txt"
TF_fname = "~/cwd/data/network/TF.txt"
V_fname = "~/cwd/data/network/V.txt"
P_eval_noknownsite = read.vector(P_eval_noknownsite_fname)
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
V = read.vector(V_fname)
PT = c(KP,TF)
nP = length(KP)
nV = length(V)
nT = length(TF)
nO = nV-length(PT)

rundir = getwd()

for (KP_edge_fname in KP_edge_fnames) {
    cat(KP_edge_fname, "\n")
    setwd(dirname(KP_edge_fname))
    KP_edge_fname = basename(KP_edge_fname)
    KP_edges = read.table(KP_edge_fname, header=T, sep="\t", quote="", check.names=F)
    P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
    # remove diagonals
    P_eval = P_eval[as.character(P_eval$Source) != as.character(P_eval$Target),]
    # sanity check
    stopifnot(all(as.character(KP_edges$KP) == as.character(P_eval$Source)))
    stopifnot(all(KP_edges$Target == P_eval$Target))
    # hold all data in eval sets
    P_eval$marker = KP_edges$marker
    top_marked = KP_edges$q < .05
    
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
    venndata$literature = venndata$biogrid | venndata$fasolo | venndata$parca | venndata$fiedler
    venndata$curated = venndata$yeastkid | venndata$ptmod
    venndata$known = venndata$literature | venndata$curated
    venndata$invitro = venndata$known | venndata$ptacek
    venndata$with_site = venndata$invitro & !P_eval$Target%in%P_eval_noknownsite
    
    p.selection = function(col_select, row_idx) {
        drawn = top_marked[row_idx]
        if(missing(row_idx)) row_idx = rep(TRUE, nrow(venndata))
        evaldata = venndata[row_idx, col_select]
        q = sum(evaldata[drawn])  # number of correct inferences
        m = sum(evaldata) # number of true edges
        n = length(evaldata) - m  # number of potential edges (not in true edge set)
        list(value=q, p=phyper(q, m, n, sum(drawn), lower.tail=FALSE))
    }
    
    KP2KP.idx = P_eval$Target%in%KP
    KP2TF.idx = P_eval$Target%in%TF
    stopifnot(all(KP2KP.idx|KP2TF.idx))
    
    ### Plotting
    KP_color = "#bd61b6"
    TF_color = "#75b42f"
    
    plot_square_euler = function(selection) {
        eulerdata.KP = data.frame(venndata[KP2KP.idx,c("potential", selection)], inferred=top_marked[KP2KP.idx])
        eulerdata.TF = data.frame(venndata[KP2TF.idx,c("potential", selection)], inferred=top_marked[KP2TF.idx])
        
        euler.plt.bg = rbind(data.frame(substrate="KP", inferred =-sum(eulerdata.KP$inferred), potential = sum(!eulerdata.KP$inferred)),
                             data.frame(substrate="TF", inferred =-sum(eulerdata.TF$inferred), potential = sum(!eulerdata.TF$inferred)))
        euler.plt.fg = rbind(data.frame(substrate="KP", TP       =-sum( eulerdata.KP$inferred & eulerdata.KP[selection]),
                                        FN       = sum(!eulerdata.KP$inferred & eulerdata.KP[selection])),
                             data.frame(substrate="TF", TP       =-sum( eulerdata.TF$inferred & eulerdata.TF[selection]),
                                        FN       = sum(!eulerdata.TF$inferred & eulerdata.TF[selection])))
        euler.plt.bg = melt(euler.plt.bg, id.vars="substrate")
        euler.plt.fg = melt(euler.plt.fg, id.vars="substrate")
        euler.plt.bg$count = abs(euler.plt.bg$value)
        # include the inferred in the actual count of potentials
        euler.plt.bg$count[euler.plt.bg$variable == "potential" & euler.plt.bg$substrate == "KP"] = sum(euler.plt.bg$count[euler.plt.bg$substrate == "KP"])
        euler.plt.bg$count[euler.plt.bg$variable == "potential" & euler.plt.bg$substrate == "TF"] = sum(euler.plt.bg$count[euler.plt.bg$substrate == "TF"])
        
        
        pvalues = c(p.selection(selection, KP2KP.idx)$p,
                    p.selection(selection, KP2TF.idx)$p)
        
        
        bw = .8
        hn = .05  # horizontal nudge
        fgscl = 6  # how wide foreground is compared to background
        
        # version where we nudge fgs down 1/4
        fgdg = bw/4
        
        ggplot(mapping=aes(x=substrate, y=value)) + 
            geom_bar(data=euler.plt.bg, mapping=aes(fill=variable), stat="identity", position="identity", width=bw, fill="transparent", color="black", size=2) +
            geom_bar(data=euler.plt.bg, mapping=aes(fill=variable), stat="identity", position="identity", width=bw) +
            geom_bar(data=euler.plt.fg, mapping=aes(fill=substrate, y=value*fgscl), stat="identity", position=position_nudge(x=-fgdg), width=bw/fgscl, alpha=.5) +
            scale_fill_manual(values=list(KP=KP_color, TF=TF_color, inferred="lightgray", potential="white")) +
            annotate("text", x=c(1,2,1,2)+bw/2-.01, vjust=1,   y=euler.plt.bg$value, hjust=c(-hn,-hn,1+hn,1+hn), label=euler.plt.bg$variable, fontface="italic", alpha=c(.65,.65,1,1)) +
            annotate("text", x=c(1,2,1,2)+bw/2-.01, vjust=2.2, y=euler.plt.bg$value, hjust=c(-hn,-hn,1+hn,1+hn), label=euler.plt.bg$count, fontface="bold", alpha=c(.65,.65,1,1)) +
            annotate("text", x=c(1,2)-fgdg+bw/fgscl, vjust=0, hjust=-hn, y=0, label="known", fontface="italic", color=c(KP_color, TF_color)) +
            annotate("text", x=c(1,2,1,2)-fgdg, vjust=0.5, y=euler.plt.fg$value*fgscl, hjust=c(1+hn,1+hn,-hn,-hn), label=abs(euler.plt.fg$value), fontface="bold", color=rep(c(KP_color, TF_color),2)) +
            annotate("text", x=c(1,2)-fgdg+bw/fgscl, vjust=0, hjust=-hn, y=10000, label=sprintf("p=%.3g",pvalues), color=c(KP_color, TF_color)) +
            coord_flip() + theme_void() + 
            theme(legend.position="none", axis.text.y=element_text(),
                  plot.margin=margin(t=0,b=0,l=0,r=0,unit="pt")) +
            scale_x_discrete(labels=parse(text=c(tex("KP\\rightarrow KP"),tex("KP\\rightarrow TF"))))
        
    }
    
    selections = c("known","literature","invitro","with_site")
    for(selection in selections) {
        ggsave(paste0("square_euler_", selection, ".pdf"), plot=plot_square_euler(selection), width=6.5, height=1.25)
    }
    
    setwd(rundir)  # go back to so relative dirs for other files still work
}


