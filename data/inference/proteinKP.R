#!/usr/bin/env Rscript

# filter KP edges for protein KP. same procedure as in ~/cwd/src/square_euler.R

library(data.table)

KP_edges = fread("~/cwd/data/inference/74/KP_edges.tsv", sep="\t", header=T)
protein_KPTF_table = read.table("~/cwd/data/network/node_attributes_full.txt", sep="\t", header=T)
protein_KP = as.character(protein_KPTF_table$ORF[protein_KPTF_table$protein_type0=="KP" & !is.na(protein_KPTF_table$protein_type0)])
protein_TF = as.character(protein_KPTF_table$ORF[protein_KPTF_table$protein_type0=="TF" & !is.na(protein_KPTF_table$protein_type0)])
protein_KPTF = c(protein_KP, protein_TF)

# filter
KP_edges = KP_edges[as.character(KP_edges$Target) %in% protein_KPTF,]
KP_edges = KP_edges[as.character(KP_edges$KP) %in% protein_KPTF,]

fwrite(KP_edges, "~/cwd/data/inference/74/KP_edges_protein.tsv", sep="\t")



