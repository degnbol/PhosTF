
# functions
flatten = function(x) as.vector(as.matrix(x))

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
# fall back on perturbation data to find expression mode
perturbation = as.matrix(read.table("../perturbation/logFC_inner.tsv", sep="\t", row.names=1, header=T, quote=""))

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
# remove the pointless data
edges = all_edges[all_edges$Pval < .5,]
edges = edges[edges$TF %in% TFs$TF,]

# create unsigned edge scores
edges$unsigned = qnorm(edges$Pval, lower.tail=F) / qnorm(0.001, lower.tail=F)
edges$unsigned[edges$unsigned == Inf] = max(edges$unsigned[edges$unsigned != Inf])
hist(edges$unsigned, breaks=30)
median(edges$unsigned)
sum(edges$unsigned > .75)
mean(edges$unsigned[edges$unsigned > .75])

# assign each TF as either activator or repressor sign by using the mode of regulation most supported in the data
# also try where it is only the unknown edges that are assigned sign
dominant_sign = function(edges) {
    signs = rep(NA, nrow(edges))
    for (TF in levels(edges$TF)) {
        activation = edges$Mode[edges$TF == TF] == "activator" & !is.na(edges$Mode[edges$TF == TF])
        inhibition = edges$Mode[edges$TF == TF] == "inhibitor" & !is.na(edges$Mode[edges$TF == TF])
        sign = activation - inhibition; sign[is.na(sign)] = 0
        # weigh by edge score
        # sign_sum = sum(sign * edges$unsigned[edges$TF == TF])
        sign_sum = sum(sign * (1-edges$Pval[edges$TF == TF]))
        if (sign_sum > 0) {signs[edges$TF == TF] = "+"}
        if (sign_sum < 0) {signs[edges$TF == TF] = "-"}
    }
    signs
}

edges$sign = dominant_sign(edges)

sum(edges$sign == "+", na.rm=T)
sum(edges$sign == "-", na.rm=T)

# remove all edges that are definitely not supported in the data, 
# we use 0.05 since it is the threshold that allows for the most edges, 
# so we have ~140000 potential edges, which is still plenty
sum(edges$Pval < .05)
edges = edges[edges$Pval < .05,]

agree = (edges$sign == "+" & edges$Mode == "activator") | (edges$sign == "-" & edges$Mode == "inhibitor")
sum(agree, na.rm=T)
sum(!agree, na.rm=T)
# maybe this is attributed to cofactors and other regulators that flip the sign of regulation

# most have been given a sign now which is great
sum(is.na(edges$sign))
sum(!is.na(edges$sign))

# provide sign to the rest from perturbation data
for(i in 1:nrow(edges)) {if (is.na(edges$sign[i])) {
    logFC = perturbation[rownames(perturbation)==edges$Target[i], colnames(perturbation)==edges$TF[i]]
    if (length(logFC) == 1) {
        if (logFC > 0) {edges$Mode[i] = "inhibitor"}
        if (logFC < 0) {edges$Mode[i] = "activator"}
    }
    else if (length(logFC) > 1) {cat("warning\n")}
}}

edges$sign[is.na(edges$sign)] = dominant_sign(edges)[is.na(edges$sign)]
# still around 4000 that is unsigned
sum(is.na(edges$sign))
sum(!is.na(edges$sign))
unique(edges$TF[is.na(edges$sign)])

hist(edges$unsigned[is.na(edges$sign)], breaks=30)

write.table(edges, file="edges.tsv", sep="\t", row.names=F, quote=F)

