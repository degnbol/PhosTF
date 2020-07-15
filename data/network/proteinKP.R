#!/usr/bin/env Rscript

library(data.table)


setwd("~/cwd/data/network")


TFs = fread("TF_mode.tsv", sep="\t")
KPs = fread("KP.tsv", sep="\t")
Vs  = fread("V.txt", header=F)$V1
measured = fread("perturbation_measured.txt", header=F)$V1
biological_DT = fread("node_attributes_full.txt", sep="\t")
TF_edges = fread("TF_priors/TF_edges.tsv", sep="\t")


# TF set is not changed
stopifnot(length(intersect(biological_DT[protein_type0 == "TF",.(ORF)]$ORF, TFs$TF)) == nrow(TFs))
stopifnot(all(TFs[TF%in%biological_DT[protein_type2=="TFA",.(ORF)]$ORF,.(Mode)]$Mode == "activator"))
stopifnot(all(TFs[TF%in%biological_DT[protein_type2=="TFR",.(ORF)]$ORF,.(Mode)]$Mode == "repressor"))

# more KPs has been classified by Chris
KPs[KP%in%biological_DT[protein_type == "PK",.(ORF)]$ORF]$Mode
KPs[KP%in%biological_DT[protein_type == "PP",.(ORF)]$ORF]$Mode

proteinKP = biological_DT[protein_type0=="KP",.(ORF, protein_type)]
fwrite(proteinKP, "KP_protein.tsv", sep="\t")

# has to be a target for regulation to be possible
Os = measured[measured%in%TF_edges$Target & !measured%in%c(proteinKP$ORF,TFs$TF)]
proteinV = c(proteinKP$ORF, TFs$TF, Os)
# there is nothing changed in the set, only small changes to the ordering
length(intersect(Vs, proteinV)) == length(proteinV)

write(proteinV, "V_protein.txt")



