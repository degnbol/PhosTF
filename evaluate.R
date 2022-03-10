#!/usr/bin/env Rscript
library(data.table)
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(verification)) # roc auc with p-value
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
# combine pvals using fishers method https://en.wikipedia.org/wiki/Fisher%27s_method 
# chisq = -2 sum(ln(p-values)); pval = 1-pchisq(chisq, df=2length(p-values))
fisher.method.log = function(pvals) {
    df = 2*length(pvals)
    pchisq(-2*sum(log(pvals),na.rm=TRUE),df,lower.tail=FALSE,log.p=TRUE)
}

# convert tex code to expression object (use \\ for \ and combine expressions with list() instead of c())
tex = function(x) unname(TeX(paste0("$",x)))


# WP_fname = "~/cwd/data/inference/02/WP_infer.mat"
WP_fnames = commandArgs(trailingOnly=T)
P_fname = "~/cwd/data/evaluation/P_eval.tsv"
KP_targets_noknownsite = read.vector("~/cwd/data/network/KP_targets_noknownsite.txt")
KP = fread("~/cwd/data/network/KP_protein.tsv")$ORF
TF = read.vector("~/cwd/data/network/TF.txt")
V = read.vector("~/cwd/data/network/V_protein.txt")
PT = c(KP,TF)
nP = length(KP)
nV = length(V)
nT = length(TF)
nO = nV-length(PT)


## settings
# assuming we are masking KP edges based on known phos site from BioGrid
masking_KP = TRUE


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
    P_eval = read.table(P_fname, header=T, sep="\t", quote="", check.names=F)
    # swap rownames and colnames columns so the order will be source then target
    KP_edges = melt_matrix(WP)[,c(2,1,3)]
    colnames(KP_edges) = c("P",  "Target", "marker")
    # sanity check
    stopifnot(all(as.character(KP_edges$P) == as.character(P_eval$Source)))
    stopifnot(all(KP_edges$Target == P_eval$Target))
    
    # hold all data in eval sets
    P_eval$marker = KP_edges$marker
    
    # remove diagonals
    P_eval = P_eval[as.character(P_eval$Source) != as.character(P_eval$Target),]
    
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
    
    # not evaluating on uninferrable edges
    if(masking_KP) {
        venndata[P_eval$Target%in%KP_targets_noknownsite,] = FALSE
    }
    
    # write KP targets that are in the evaluation dataset. This was used before to make KP_targets_noknownsite, not sure if we still use it now.
    # write.table(PT[PT%in%P_eval$Target[venndata$invitro]], "~/cwd/data/evaluation/KP_targets.txt", row.names=F, col.names=F, quote=F)
    
    # index indicating which edge is inferred most strongly to most weakly
    marker_order = function(k, index) {
        if(missing(index)) marker = P_eval$marker
        else marker = P_eval$marker[index]
        order(abs(marker), decreasing=T)[1:k]
    }
    
    # logical version of marker_order
    top_marked = function(k, index) {
        if(!missing(index)) dim=sum(index)
        else dim=nrow(P_eval)
        unwhich(marker_order(k, index), dim=dim)
    }
    
    p.selection = function(k, col_select, row_idx) {
        drawn = top_marked(k, row_idx)
        if(missing(row_idx)) row_idx = rep(TRUE, nrow(venndata))
        evaldata = venndata[row_idx, col_select]
        q = sum(evaldata[drawn])  # number of correct inferences
        m = sum(evaldata) # number of true edges
        n = sum(venndata[row_idx, "potential"]) - m  # number of black balls to pick/"false" edges/potential edges not in set of known
        list(value=q, p=phyper(q, m, n, sum(drawn), lower.tail=FALSE))
    }
    
    quantiles = seq(0.01, 0.25, 0.0025)
    quantiles_plot = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.25)
    
    p.selection.q = function(quantiles, col_select, row_idx) {
        out = data.frame()
        if (missing(row_idx)) n = nrow(venndata)
        else n = sum(row_idx)
        for (k in quantiles*n) {
            out = rbind(out, p.selection(k, col_select, row_idx))
        }
        data.frame(test=paste(col_select, quantiles), selection=col_select, quantile=quantiles, out)
    }
    
    KP2KP.idx = P_eval$Target%in%KP
    KP2TF.idx = P_eval$Target%in%TF
    stopifnot(all(KP2KP.idx|KP2TF.idx))
    
    eval.table = data.frame()
    selections = c("known","literature","invitro")
    for (selection in selections) {
        selection.table = rbind(data.frame(p.selection.q(quantiles, selection, KP2KP.idx), substrate="KP"),
                                data.frame(p.selection.q(quantiles, selection, KP2TF.idx), substrate="TF"))
        eval.table = rbind(eval.table, selection.table)
    }
    
    # only select few quantiles in the eval table
    eval.table.few = eval.table[eval.table$quantile%in%quantiles_plot,]
    eval.table.few = eval.table.few[order(eval.table.few$p),]  # sort by p-value
    write.table(eval.table.few, "evaluation.tsv", sep="\t", quote=F, row.names=F)
    
    # best selection
    fisher.ps = aggregate(p ~ selection, data=eval.table, fisher.method.log)
    best.selection = fisher.ps$selection[which.min(fisher.ps$p)]
    
    write.table(-fisher.ps$p[best.selection], "score.txt", quote=F, row.names=F, col.names=F)
    
    ### Plotting
    KP_color = "#bd61b6"
    TF_color = "#75b42f"
    
    plot_square_euler = function(selection) {
        best.quantile = .15
        inferred.KP = top_marked(best.quantile * sum(KP2KP.idx), KP2KP.idx)
        inferred.TF = top_marked(best.quantile * sum(KP2TF.idx), KP2TF.idx)
        eulerdata.KP = data.frame(venndata[KP2KP.idx,c("potential", selection)], inferred=inferred.KP)
        eulerdata.TF = data.frame(venndata[KP2TF.idx,c("potential", selection)], inferred=inferred.TF)
        
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
            coord_flip() + theme_void() + 
            theme(legend.position="none", axis.text.y=element_text(),
                  plot.margin=margin(t=0,b=0,l=0,r=0,unit="pt")) +
            scale_x_discrete(labels=parse(text=c(tex("KP\\rightarrow KP"),tex("KP\\rightarrow TF"))))
        
    }
    
    
    plot_logp = function(plt.p) {
        
        second_axis = dup_axis(name="substrates/KP", breaks=quantiles_plot, labels=round(quantiles_plot*nrow(venndata)/length(KP),1))
        
        ggplot(data=plt.p, aes(x=quantile, y=-log10(p), color=substrate)) +
            geom_line() +
            scale_y_continuous(breaks=seq(0,30,3), expand=c(0,0), limits=c(0, 30)) +
            scale_x_continuous(breaks=quantiles_plot, labels=gsub("0\\.","\\.",quantiles_plot), sec.axis=second_axis, expand=c(0,0), limits=c(min(quantiles_plot),max(quantiles_plot))) +
            theme_linedraw() +
            geom_hline(yintercept=-log10(0.05), linetype="dashed") + 
            annotate("text", x=.31, y=2.5, label="p=0.05", size=3) +
            ylab(expression("-"*log[10]*" p")) +
            theme(panel.grid.major=element_line(colour="gray"), 
                  panel.grid.minor=element_line(colour="lightgray"),
                  plot.margin=margin(t=0,b=0,l=0,r=0,unit="pt"),
                  legend.title=element_blank()) +
            scale_color_manual(values=c(KP_color, TF_color), labels=list(tex("KP\\rightarrow KP"),tex("KP\\rightarrow TF")))
    }
    
    for(selection in selections) {
        plt.p = eval.table[eval.table$selection==selection,]
        ggsave(paste0("square_euler_",selection,".pdf"), plot=plot_square_euler(selection), width=6.5, height=1.25)
        ggsave(paste0("hyperp_",selection,".pdf"), plot=plot_logp(plt.p), width=6.5, height=2.5)
    }
    
    setwd(rundir)  # go back to so relative dirs for other files still work
}


