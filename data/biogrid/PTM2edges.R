#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data/biogrid"))
# REQUIRES: raw/BIOGRID-PTM.ptmtab.txt and raw/BIOGRID-PTM-RELATIONSHIPS.ptmrel.txt

# read
targets = fread("raw/BIOGRID-PTM.ptmtab.txt", sep="\t", header=TRUE, quote="")
sources = fread("raw/BIOGRID-PTM-RELATIONSHIPS.ptmrel.txt", sep="\t", header=TRUE, quote="")
setnames(sources, "#PTM ID", "PTM.ID")
setnames(targets, "#PTM ID", "PTM.ID")
mapping = match(sources$PTM.ID, targets$PTM.ID)
# combine
names(targets) = paste("Target", names(targets), sep=".")
edges = cbind(sources, targets[mapping,])
# remove NA
edges_nona = na.omit(edges)
cat((1- nrow(edges_nona) / nrow(edges)) * 100, "percent of edges were removed\n")
# sanity check
stopifnot(all(edges_nona$PTM.ID==edges_nona$Target.PTM.ID))
# remove duplicate PTM ID column
edges_nona[,Target.PTM.ID:=NULL]
names(edges_nona) = gsub("\\.", " ", names(edges_nona))
fwrite(edges_nona, file="regulatory_edges.tsv", sep="\t", quote=F)
edges_core = edges_nona[`Target Post Translational Modification`=="Phosphorylation", .(P=`Systematic Name`, Target=`Target Systematic Name`, Relationship)]
fwrite(edges_core, file="P_edges.tsv", sep="\t", quote=F)

