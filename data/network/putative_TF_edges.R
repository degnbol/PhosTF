
library(data.table)
library(ggplot2)

# functions

# combine pvals using fishers method https://en.wikipedia.org/wiki/Fisher%27s_method 
# chisq = -2 sum(ln(p-values)); pval = 1-pchisq(chisq, df=2length(p-values))
fisher.method = function(pvals) {
    df = 2*length(pvals)
    pchisq(-2*sum(log(pvals),na.rm=TRUE),df,lower.tail=FALSE)
}


setwd("~/cwd/data/network")

# read
balaji_2006 = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", header=T)[,1:2]
# collected from harbison, lee 2002 (merged into harbison), and horak. so only use as upper bound for pval (pval<0.001)
balaji_2008 = read.table("../processed/balaji_2008/TF_edges.tsv", sep="\t", header=T)
harbison_YPD = read.table("../processed/harbison_2004/TF_edges_YPD.tsv", sep="\t", header=T)
harbison_conds = read.table("../processed/harbison_2004/TF_edges_conds.tsv", sep="\t", header=T)
horak = read.table("../processed/horak_2002/TF_edges.tsv", sep="\t", header=T)[,c(1,2,5)]; colnames(horak)[3] = "Pval"
workman = read.table("../processed/workman_2006/TF_edges.tsv", sep="\t", header=T)[,c(1,2,4)]
yeastract_binding = read.table("../processed/yeastract/binding_ORF.tsv", sep="\t", col.names=c("TF", "Target"))
# expression data used to provide sign to edges, not to indicate presence of edges
yeastract_expression = read.table("../processed/yeastract/expression_ORF.tsv", sep="\t", col.names=c("TF", "Target", "Mode"), na.strings="ambiguous")
yeastract_expression = na.omit(yeastract_expression)
STRING = read.table("../processed/STRING/modes.tsv", sep="\t", header=T); colnames(STRING)[1] = "TF"
# combine expression evidence, trust yeastract more than STRING since yeastract only has non-conflicting evidence, 
# while STRING is a combined score and in cases where there is conflict it uses the highest scoring (activator or repressor)
expression_evidence = rbind(yeastract_expression, STRING[!(paste(STRING$TF,STRING$Target) %in% paste(yeastract_expression$TF,yeastract_expression$Target)),])

all_edges = rbind(harbison_YPD, harbison_conds, horak, workman)
# pval=0 doesn't make sense, change it to the smallest nonzero one
all_edges$Pval[all_edges$Pval==0] = min(all_edges$Pval[all_edges$Pval!=0])
# min(all_edges$Pval) == 1.11e-16

# add edges from reviews with p-val thresholds
balaji = unique(rbind(balaji_2006, balaji_2008))
# remove edges from balaji that are already found in original experiment
balaji = balaji[!(paste(balaji$TF, balaji$Target) %in% paste(all_edges$TF, all_edges$Target)),]
balaji$Pval = 0.001
all_edges = rbind(all_edges, balaji)
# remove edges yeastract that are already found
yeastract_binding = yeastract_binding[!(paste(yeastract_binding$TF, yeastract_binding$Target) %in% paste(all_edges$TF, all_edges$Target)),]
yeastract_binding$Pval = 0.04999 # I couldn't find a pval threshold so I will use the worst significance level they could have chosen
all_edges = rbind(all_edges, yeastract_binding)

# apply fisher's method to combine pvals for each edge
all_edges = data.table(all_edges)[,list(Pval=fisher.method(Pval)),by=list(TF,Target)]

# add regulation mode supported in the data
all_edges = merge(all_edges, expression_evidence, all.x=T)

plot(density(all_edges$Pval))
# remove all edges not supported in the data, we use 0.05 since it is the threshold that allows for the most edges.
# we get ~80000 potential edges, which is still plenty
sum(all_edges$Pval <= .05)
edges = all_edges[all_edges$Pval <= .05,]

write.table(edges, "TF_edges_putative.tsv", sep="\t", row.names=F, quote=F, na="")



