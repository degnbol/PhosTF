
# purpose of this is to assign mode to the TFs that did not have GO terms that could do it.
# these modes are therefore estimates, but we need to assign something, so that their outgoing edges can be signed
# mode inferred from datasets with mode annotation for edges (yeastract) combined with the belief in the edge (based on edge p-val)

setwd("/Users/christian/cwd/data/nodes")

TFs = read.table("../nodes/TF.tsv", sep="\t", header=T, stringsAsFactors=F)
V = as.vector(as.matrix(read.table("../nodes/V.txt")))
edges = read.table("../TF_edges/edges.tsv", sep="\t", header=T, quote="", stringsAsFactors=F)
# fall back on perturbation data to find expression mode
perturbation = as.matrix(read.table("../perturbation/logFC_inner.tsv", sep="\t", row.names=1, header=T, quote="", stringsAsFactors=F))

# assign each TF as either activator or repressor sign by using the mode of regulation most supported in the data
dominant_sign = function(edges) {
    # pval=0 gives Inf so
    edges$Pval[edges$Pval == 0] = min(edges$Pval[edges$Pval != 0])
    signs = rep(NA, nrow(edges))
    for (TF in unique(edges$TF)) {
        activation = edges$Mode[edges$TF == TF] == "activator"
        inhibition = edges$Mode[edges$TF == TF] == "inhibitor"
        chisq_activation = -2*sum(log(edges$Pval[activation]))
        chisq_inhibition = -2*sum(log(edges$Pval[inhibition]))
        pval_activation = pchisq(chisq_activation, df=2*sum(activation), lower.tail=F, log.p=T)
        pval_inhibition = pchisq(chisq_inhibition, df=2*sum(inhibition), lower.tail=F, log.p=T)
        
        if (pval_activation < pval_inhibition) {signs[edges$TF == TF] = "+"}
        if (pval_activation > pval_inhibition) {signs[edges$TF == TF] = "-"}
    }
    signs
}

edges$sign = dominant_sign(edges)

sum(edges$sign == "+", na.rm=T)
sum(edges$sign == "-", na.rm=T)


agree = (edges$sign == "+" & edges$Mode == "activator") | (edges$sign == "-" & edges$Mode == "inhibitor")
sum(agree)
sum(!agree & edges$Mode != "")
# maybe this is attributed to cofactors and other regulators that flip the sign of regulation

# most have been given a sign now which is great
sum(is.na(edges$sign)) / sum(!is.na(edges$sign))


# provide sign from perturbation data where a TF is KOed
for(i in 1:nrow(edges)) {if (is.na(edges$sign[i])) {
    logFC = perturbation[rownames(perturbation)==edges$Target[i], colnames(perturbation)==edges$TF[i]]
    if (length(logFC) == 1) {
        if (logFC > 0) {edges$Mode[i] = "inhibitor"}
        if (logFC < 0) {edges$Mode[i] = "activator"}
    }
    else if (length(logFC) > 1) {cat("warning\n")}
}}

edges$sign[is.na(edges$sign)] = dominant_sign(edges)[is.na(edges$sign)]
unresolved = unique(edges$TF[is.na(edges$sign)])
length(unresolved) # 19 TFs that are unsigned (15, 4 are known from GO terms)

TF_signs = edges$sign[match(TFs$TF, edges$TF)]

# get the sign of expression correlation between an unresolved TF and its targets
for (TF in unresolved) {if (TF %in% rownames(perturbation)) {
    targets = edges$Target[edges$TF == TF]
    targets = targets[targets %in% rownames(perturbation)]
    corr = as.vector(cor(rep(perturbation[TF,], length(targets)), matrix(perturbation[targets,], byrow=T)))
    if (corr < 0) {TF_signs[TFs$TF == TF] = "-"}
    if (corr > 0) {TF_signs[TFs$TF == TF] = "+"}
}}


TFs$Mode[TF_signs == "+" & !is.na(TF_signs) & TFs$Mode == ""] = "activator"
TFs$Mode[TF_signs == "-" & !is.na(TF_signs) & TFs$Mode == ""] = "repressor"


# manually looking up the only two remaining TFs: MATA1 and YER108C
# SGD says MATA1 represses genes, and YER108C is a point mutated copy of YER109C (FLO8) which is an activator
TFs$Mode[TFs$TF == "MATA1"] = "repressor"
TFs$Mode[TFs$TF == "YER108C"] = "activator"

write.table(TFs, "TF_mode.tsv", sep="\t", quote=F, row.names=F)
