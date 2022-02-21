#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data/biogrid"))
# REQUIRES: yeast_PTM.tsv and yeast_PTM_relationships.tsv

# read
targets = fread("yeast_PTM.tsv", sep="\t", header=TRUE, quote="")
sources = fread("yeast_PTM_relationships.tsv", sep="\t", header=TRUE, quote="")
setnames(sources, "#PTM ID", "PTM.ID")
setnames(targets, "#PTM ID", "PTM.ID")
mapping = match(sources$PTM.ID, targets$PTM.ID)
# combine
names(targets) = paste("Target", names(targets), sep=".")
edges = cbind(sources, targets[mapping,])
# remove NA
edges_nona = na.omit(edges)
cat((1- nrow(edges_nona) / nrow(edges)) * 100, "percent of edges had NA so were removed\n")
# sanity check
stopifnot(all(edges_nona$PTM.ID==edges_nona$Target.PTM.ID))
# remove duplicate PTM ID column
edges_nona[,Target.PTM.ID:=NULL]
names(edges_nona) = gsub("\\.", " ", names(edges_nona))
fwrite(edges_nona, file="regulatory_edges.tsv", sep="\t", quote=F)
edges_phos = edges_nona[`Target Post Translational Modification`=="Phosphorylation", .(P=`Systematic Name`, Target=`Target Systematic Name`, Relationship)]
edges_phos = unique(edges_phos)
fwrite(edges_phos, file="P_edges.tsv", sep="\t", quote=F)

