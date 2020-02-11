
library(data.table)
library(ggplot2)
library(harmonicmeanp)

# putative edge refers to an edge with significant p-value (<=0.05), where it is putative if the source of the edge is actually a TF.
# this is settled in the ../nodes.R script which should be run next.


# functions

# combine pvals using fishers method https://en.wikipedia.org/wiki/Fisher%27s_method 
# chisq = -2 sum(ln(p-values)); pval = 1-pchisq(chisq, df=2length(p-values))
fisher.method = function(pvals) {
    df = 2*length(pvals)
    pchisq(-2*sum(log(pvals),na.rm=TRUE),df,lower.tail=FALSE)
}


setwd("~/cwd/data/network/TF_priors")

# read
balaji_2006 = read.table("~/cwd/data/processed/balaji_2006/TF_edges.tsv", sep="\t", header=T)[,1:2]
# collected from harbison, lee 2002 (merged into harbison), and horak. so only use as upper bound for pval (pval<0.001)
balaji_2008 = read.table("~/cwd/data/processed/balaji_2008/TF_edges.tsv", sep="\t", header=T)
harbison_YPD = read.table("~/cwd/data/processed/harbison_2004/TF_edges_YPD.tsv", sep="\t", header=T)
harbison_conds = read.table("~/cwd/data/processed/harbison_2004/TF_edges_conds.tsv", sep="\t", header=T)
horak = read.table("~/cwd/data/processed/horak_2002/TF_edges.tsv", sep="\t", header=T)[,c(1,2,5)]; colnames(horak)[3] = "Pval"
workman = read.table("~/cwd/data/processed/workman_2006/TF_edges.tsv", sep="\t", header=T)[,c(1,2,4)]
yeastract_binding = read.table("~/cwd/data/processed/yeastract/binding_ORF.tsv", sep="\t", col.names=c("TF", "Target"))
# expression data used to provide sign to edges, not to indicate presence of edges
yeastract_expression = read.table("~/cwd/data/processed/yeastract/expression_ORF.tsv", sep="\t", col.names=c("TF", "Target", "Mode"), na.strings="ambiguous")
yeastract_expression = na.omit(yeastract_expression)
STRING = read.table("~/cwd/data/processed/STRING/modes.tsv", sep="\t", header=T); colnames(STRING)[1] = "TF"
# combine expression evidence, trust yeastract more than STRING since yeastract only has non-conflicting evidence, 
# while STRING is a combined score and in cases where there is conflict it uses the highest scoring (activator or repressor)
expression_evidence = rbind(yeastract_expression, STRING[!(paste(STRING$TF,STRING$Target) %in% paste(yeastract_expression$TF,yeastract_expression$Target)),])

all_edges = rbind(harbison_YPD, harbison_conds, horak, workman)
# pval=0 doesn't make sense, change it to the smallest nonzero one
all_edges$Pval[all_edges$Pval==0] = min(all_edges$Pval[all_edges$Pval!=0])
# min(all_edges$Pval) == 1.11e-16
# apply fisher's method to combine pvals assumed independent for each edge
all_edges = data.table(all_edges)[,list(Pval=fisher.method(Pval)),by=list(TF,Target)]

balaji = unique(rbind(balaji_2006, balaji_2008))
balaji$Pval = 0.001
# I couldn't find a pval threshold so I will use the worst significance level that will still be included in the end
yeastract_binding$Pval = 0.04999

all_edges = rbind(all_edges, balaji, yeastract_binding)
# use data that is not independent to set a minimum p-value for the data
all_edges = all_edges[,list(Pval=min(Pval)),by=list(TF,Target)]

# add regulation mode supported in the data
all_edges = merge(all_edges, expression_evidence, all.x=T)

hist(all_edges$Pval, breaks=100)

library(fdrtool)

# remove all edges not supported in the data, we use 0.05 since it is the threshold that allows for the most edges.
# we get ~100000 potential edges, which is still plenty
sum(all_edges$Pval <= .05)
sum(p.adjust(all_edges$Pval, "BH") <= .05)
fdr = fdrtool(all_edges$Pval, statistic="pvalue", plot=FALSE)
sum(fdr$qval < .2) / 231
write.table(all_edges[fdr$qval <= .2,], "TF_edges_putative_FDR.tsv", sep="\t", row.names=F, quote=F, na="")
write.table(all_edges[all_edges$Pval <= .05,], "TF_edges_putative.tsv", sep="\t", row.names=F, quote=F, na="")



