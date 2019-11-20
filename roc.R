# packages
library(ggplot2)
library(plotROC)
library(dplyr)
library(extrafont)

# functions
read_vec = function(fname) {as.vector(as.matrix(read.table(fname)))}
get_pos = function(v) {out = +v; out[out <= 0] = 0; out}
get_neg = function(v) {out = -v; out[out <= 0] = 0; out}

# filenames
setwd("/Users/christian/GoogleDrev/PKTFX/data/dream/archive/LB01_100")
true_fname = "WP.mat"
marker_fname = "WP_infer.mat"

# read
true_fnames = list.files(pattern=true_fname, recursive=T)
marker_fnames = list.files(pattern=marker_fname, recursive=T)
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
    geom_roc(n.cuts=0, aes(fill=sample), linealpha=.2) + 
    geom_roc(n.cuts=0) + 
    style_roc(theme=theme_bw, guide=F, xlab="FPR", ylab="TPR") +
    scale_color_manual(labels=labels, values=c("red", "green", "blue")) +
    theme(legend.text=element_text(family="Menlo"), legend.title=element_blank(), 
          legend.position=c(1,0), legend.justification=c(1,0), legend.box.margin=margin(1, 1, 1, 1))

plot


