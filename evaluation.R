#!/usr/bin/env Rscript

# packages
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


WP_fnames = commandArgs(trailingOnly=T)
# WP_fnames = "~/cwd/data/inference/02/WP_infer.mat"
P_fname = "~/cwd/data/evaluation/P_eval.tsv"
KP_fname = "~/cwd/data/network/KP.txt"
TF_fname = "~/cwd/data/network/TF.txt"
V_fname = "~/cwd/data/network/V.txt"
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
V = read.vector(V_fname)
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
    
    evaluate_cor = function(dataset, cor_names, cor_names_pos, cor_names_neg) {
        out = data.frame()
        for (name in c(cor_names, cor_names_pos, cor_names_neg)) {
            valid = !is.na(dataset[name])
            if (name%in%cor_names)   marker = abs(dataset[valid, "marker"])
            if (name%in%cor_names_pos) marker = +(dataset[valid, "marker"])
            if (name%in%cor_names_neg) marker = -(dataset[valid, "marker"])
            PCC = cor.test(dataset[valid,name], marker)[c("estimate", "p.value")]
            out = rbind(out, data.frame(test=paste(name, "cor"), value=PCC[[1]][[1]], p=PCC[[2]]))
        }
        out
    }
    
    cor_names = c("yeastkid", "reaction", "ptmod", "expression", "catalysis", "netphorest", "networkin", "networkin STRING", 
                  "networkin_biogrid", "undirected", "EMAP", "n_datasets")
    cor_names_pos = c("activation")
    cor_names_neg = c("inhibition")
    cortable = evaluate_cor(P_eval, cor_names, cor_names_pos, cor_names_neg)
    
    aucs = rbind(
        data.frame(roc.area(P_eval$n_datasets > 0, abs(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(!is.na(P_eval$parca), abs(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(!is.na(P_eval$fasolo), abs(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(P_eval$fiedler != "", abs(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(P_eval$biogrid != "", abs(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(P_eval$biogrid == "kinase", +(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(P_eval$biogrid == "phosphatase", -(P_eval$marker))[c("A","p.value")]),
        data.frame(roc.area(P_eval$goldstandard1[!is.na(P_eval$goldstandard1)], abs(P_eval$marker[!is.na(P_eval$goldstandard1)]))[c("A","p.value")]),
        data.frame(roc.area(P_eval$goldstandard2[!is.na(P_eval$goldstandard2)], abs(P_eval$marker[!is.na(P_eval$goldstandard2)]))[c("A","p.value")]),
        data.frame(roc.area(P_eval$goldstandard3[!is.na(P_eval$goldstandard3)], abs(P_eval$marker[!is.na(P_eval$goldstandard3)]))[c("A","p.value")]),
        data.frame(roc.area(P_eval$goldstandard4[!is.na(P_eval$goldstandard4)], abs(P_eval$marker[!is.na(P_eval$goldstandard4)]))[c("A","p.value")])
    )
    titles = c("any", "parca", "fasolo", "fiedler", "biogrid", "kinase", "phosphatase", "gold1", "gold2", "gold3", "gold4")
    aucs_P = cbind(test=paste(titles, "auc"), aucs)
    colnames(aucs_P)[2:3] = c("value", "p")
    
    
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
        n = length(evaldata) - m  # number of potential edges (not in true edge set)
        list(value=q, p=phyper(q, m, n, sum(drawn), lower.tail=FALSE))
    }
    
    quantiles = seq(0.01, 0.333, 0.0025)
    quantiles_plot = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.25, 0.333)
    
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
    
    known_KP.p = p.selection.q(quantiles, "known", KP2KP.idx)
    literature_KP.p = p.selection.q(quantiles, "literature", KP2KP.idx)
    invitro_KP.p = p.selection.q(quantiles, "invitro", KP2KP.idx)
    known_TF.p = p.selection.q(quantiles, "known", KP2TF.idx)
    literature_TF.p = p.selection.q(quantiles, "literature", KP2TF.idx)
    invitro_TF.p = p.selection.q(quantiles, "invitro", KP2TF.idx)
    known_KP.p$substrate = literature_KP.p$substrate = invitro_KP.p$substrate = "KP"
    known_TF.p$substrate = literature_TF.p$substrate = invitro_TF.p$substrate = "TF"
    known.p = rbind(known_KP.p, known_TF.p)
    literature.p = rbind(literature_KP.p, literature_TF.p)
    invitro.p = rbind(invitro_KP.p, invitro_TF.p)
    positive.intersections = rbind(known.p, literature.p, invitro.p)
    eval.table = rbind(cortable, aucs_P, positive.intersections[positive.intersections$quantile%in%quantiles_plot, c("test", "value", "p")])
    eval.table = eval.table[order(eval.table$p),]  # sort by p-value
    
    fisher.ps = list(known=fisher.method.log(c(known_KP.p$p, known_TF.p$p)),
                     literature=fisher.method.log(c(literature_KP.p$p, literature_TF.p$p)),
                     invitro=fisher.method.log(c(invitro_KP.p$p, invitro_TF.p$p)))
    # best selection 
    best.selection = names(which.min(fisher.ps))
    plt.p = list(known=known.p, literature=literature.p, invitro=invitro.p)[[best.selection]]
    # best quantile
    both.p = aggregate(p ~ quantile, data=plt.p, fisher.method.log)
    best.quantile = both.p[order(both.p$p),][1,"quantile"]
    
    write.table(eval.table, "evaluation.tsv", sep="\t", quote=F, row.names=F)
    write.table(-fisher.ps[[best.selection]], "score.txt", quote=F, row.names=F, col.names=F)
    
    ### Plotting
    KP_color = "#bd61b6"
    TF_color = "#75b42f"
    
    selection = as.character(best.selection)
    inferred.KP = top_marked(best.quantile * sum(KP2KP.idx), KP2KP.idx)
    inferred.TF = top_marked(best.quantile * sum(KP2TF.idx), KP2TF.idx)
    eulerdata.KP = data.frame(venndata[KP2KP.idx,c("potential", selection)], inferred=inferred.KP)
    eulerdata.TF = data.frame(venndata[KP2TF.idx,c("potential", selection)], inferred=inferred.TF)
    
    plot_eulers = function() {
        eulerr_options(padding=unit(14, units="pt"), labels=list(font=3), quantities=list(font=2))
        pdf("euler_KP.pdf", width=4, height=4)
        print(plot(euler(eulerdata.KP, shape="ellipse"), quantities=T, labels=c("potential", "known", "inferred"), adjust_labels=T, fills=list(fill=c("transparent", KP_color, "lightgray"))))
        dev.off()
        pdf("euler_TF.pdf", width=4, height=4)
        print(plot(euler(eulerdata.TF, shape="ellipse"), quantities=T, labels=c("potential", "known", "inferred"), adjust_labels=T, fills=list(fill=c("transparent", TF_color, "lightgray"))))
        dev.off()
    }
    
    plot_square_euler = function() {
        
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
        # ggplot(mapping=aes(x=substrate, y=value)) + 
        #     geom_bar(data=euler.plt.bg, mapping=aes(fill=variable), stat="identity", position="identity", width=bw, fill="transparent", color="black", size=2) +
        #     geom_bar(data=euler.plt.bg, mapping=aes(fill=variable), stat="identity", position="identity", width=bw) +
        #     geom_bar(data=euler.plt.fg, mapping=aes(fill=substrate, y=value*fgscl), stat="identity", position="identity", width=bw/fgscl, alpha=.5) +
        #     scale_fill_manual(values=list(KP=KP_color, TF=TF_color, inferred="lightgray", potential="white")) +
        #     annotate("text", x=c(1,2,1,2)+bw/2-.01, vjust=1,   y=euler.plt.bg$value, hjust=c(-hn,-hn,1+hn,1+hn), label=euler.plt.bg$variable, fontface="italic", alpha=c(.65,.65,1,1)) +
        #     annotate("text", x=c(1,2,1,2)+bw/2-.01, vjust=2.2, y=euler.plt.bg$value, hjust=c(-hn,-hn,1+hn,1+hn), label=euler.plt.bg$count, fontface="bold", alpha=c(.65,.65,1,1)) +
        #     annotate("text", x=c(1,2)+bw/4, vjust=1, hjust=0-hn, y=0, label="known", fontface="italic", color=c(KP_color, TF_color)) +
        #     annotate("text", x=c(1,2,1,2), vjust=0.5, y=euler.plt.fg$value*fgscl, hjust=c(1+hn,1+hn,-hn,-hn), label=abs(euler.plt.fg$value), fontface="bold", color=rep(c(KP_color, TF_color),2)) +
        #     coord_flip() + theme_void() + 
        #     theme(legend.position="none", axis.text.y=element_text()) +
        #     scale_x_discrete(labels=parse(text=c(tex("KP\\rightarrow KP"),tex("KP\\rightarrow TF"))))
        
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
    
    
    plot_logp = function() {

        second_axis = dup_axis(name="substrates/KP", breaks=quantiles_plot, labels=round(quantiles_plot*nrow(venndata)/length(KP),1))
        
        ggplot(data=plt.p, aes(x=quantile, y=-log10(p), color=substrate)) +
            geom_line() +
            scale_y_continuous(breaks=seq(0,27,3), expand=c(0,0), limits=c(0, 27)) +
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
    
    
    ggsave("square_euler.pdf", plot=plot_square_euler(), width=6.5, height=1.25)
    ggsave("hyperp.pdf", plot=plot_logp(), width=6.5, height=2.5)
    
    setwd(rundir)  # go back to so relative dirs for other files still work
}

