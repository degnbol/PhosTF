#!/usr/bin/env Rscript

library(data.table)

#######################################################
# count edges KP->KP, KP->TF, TF->KP, TF->TF
#######################################################

args = commandArgs(T)
WP = as.matrix(fread(args[1]))
WT = as.matrix(fread(args[2]))

nKP = ncol(WP)

nKP2KP = sum(WP[1:nKP,] != 0)
nKP2TF = sum(WP[nKP+1:nrow(WP),] != 0)
nTF2KP = sum(WT[1:nKP,] != 0)
nTF2TF = sum(WT[nKP+1:nrow(WT),] != 0)

cat(paste(c(nKP2KP,nKP2TF,nTF2KP,nTF2TF), sep="\n"), "\n")


