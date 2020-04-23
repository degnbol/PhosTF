#!/usr/bin/env Rscript

library(data.table)

#######################################################
# count edges KP->KP, KP->TF, TF->KP, TF->TF, TF->O
#######################################################

# first half is WPs second half is WTs
args = commandArgs(T)
stopifnot(length(args) %% 2 == 0)
n = length(args) / 2

for (i in 1:n) {
    WP = as.matrix(fread(args[i]))
    WT = as.matrix(fread(args[i+n]))
    
    nKP = ncol(WP)
    nTF = ncol(WT)
    
    nKP2KP = sum(WP[1:nKP,] != 0)
    nKP2TF = sum(WP[(nKP+1):nrow(WP),] != 0)
    nTF2KP = sum(WT[1:nKP,] != 0)
    nTF2TF = sum(WT[nKP+(1:nTF),] != 0)
    nTF2O  = sum(WT[(nKP+nTF+1):nrow(WT),] != 0)
    
    cat(paste(c(nKP2KP,nKP2TF,nTF2KP,nTF2TF,nTF2O)), "\n")
}

