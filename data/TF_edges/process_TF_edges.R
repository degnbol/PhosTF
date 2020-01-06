
setwd("~/cwd/data/TF_edges")

# read
TFs = read.table("../nodes/TF_mode.tsv", sep="\t", header=T)
balaji = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", header=T)
balaji_2008 = read.table("../processed/balaji_2008/TF_edges.tsv", sep="\t", header=T)
harbison = read.table("../processed/harbison_2004/TF_edges.tsv", sep="\t", header=T); colnames(harbison)[2] = "Target"
horak = read.table("../processed/horak_2002/TF_edges.tsv", sep="\t", header=T)[,c(1,2,5)]; colnames(horak)[3] = "Pval"
workman = read.table("../processed/workman_2006/TF_edges.tsv", sep="\t", header=T)[,c(1,2,4)]
yeastract_binding = read.table("../processed/yeastract/binding_ORF.tsv", sep="\t", col.names=c("TF", "Target"))
# expression data used to provide sign to edges, not to indicate presence of edges
yeastract_expression = read.table("../processed/yeastract/expression_ORF.tsv", sep="\t", col.names=c("TF", "Target", "Mode"), na.strings="ambiguous")
yeastract_expression = na.omit(yeastract_expression)
STRING = read.table("../processed/STRING/modes.tsv", sep="\t", header=T); colnames(STRING)[1] = "TF"

# first collect edges with scores
all_edges = aggregate(Pval ~ ., data=rbind(harbison, horak, workman), min)
# assign median pvals where pvals are not given for each edge
balaji$Pval = balaji_2008$Pval = median(all_edges$Pval[all_edges$Pval < 0.001])
yeastract_binding$Pval = median(all_edges$Pval[all_edges$Pval < 0.05])
# add the data that now has median pvals
all_edges = aggregate(Pval ~ ., data=rbind(all_edges, balaji, balaji_2008, yeastract_binding), min)

# add regulation mode supported in the data
all_edges = merge(all_edges, yeastract_expression, all.x=T)
all_edges = merge(all_edges, STRING, all.x=T)

plot(density(all_edges$Pval)) # peak at the far left starts at 0.05
# remove all edges that are definitely not supported in the data, 
# we use 0.05 since it is the threshold that allows for the most edges, 
# so we have ~140000 potential edges, which is still plenty
sum(all_edges$Pval < .05)
edges = all_edges[all_edges$Pval < .05,]
edges = edges[edges$TF %in% TFs$TF,]
edges = edges[as.character(edges$TF) != as.character(edges$Target),] # remove self-loops
hist(edges$Pval)

# create unsigned edge scores
edges$unsigned = qnorm(edges$Pval, lower.tail=F) / qnorm(1e-4, lower.tail=F)
edges$unsigned[edges$unsigned == Inf] = max(edges$unsigned[edges$unsigned != Inf])
hist(edges$unsigned, breaks=30)
median(edges$unsigned)
c(sum(edges$unsigned >= .5), sum(edges$unsigned >= .75), sum(edges$unsigned >= 1.))
mean(edges$unsigned[edges$unsigned > .75])

write.table(edges, file="edges.tsv", sep="\t", row.names=F, quote=F, na="")
