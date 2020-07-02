#!/usr/bin/env Rscript

library(data.table)


KP_edges = fread("~/cwd/data/inference/74/KP_edges_eval_noTN_protein.tsv", sep="\t", header=T)
perturbation = fread("~/cwd/data/perturbation/logFC_inner.csv", sep=",", header=T)
# make sure we are reading a version where first column is named ORF (data.table format) and not without row names header (data.frame format)
for(n in colnames(perturbation)) {stopifnot(is.character(n))}

perturb_names = colnames(perturbation)[2:ncol(perturbation)]

KP_edges = KP_edges[(target_type=="TF") & (inferred == T) & (KP%in%perturb_names) & (Target%in%perturb_names),]

pert_prods = function(edges) {
    out = c()
    for (i in 1:nrow(edges)) {
        out = c(out, perturbation[,edges$KP[i],with=F] * perturbation[, edges$Target[i], with=F])
    }
    return(unlist(out))
}

prods = pert_prods(KP_edges[confusion=="TP",])


sum(prods < -1)
sum(prods > +1)




