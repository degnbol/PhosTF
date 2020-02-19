#!/usr/bin/env Rscript
setwd("~/cwd/data/network/KP_priors")
load("wilcoxon.RData")
# packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

kps = commandArgs(trailingOnly=TRUE)
cat(length(kps),"\n")

for(kp in kps) {
    cat(kp, "\n")
    # all combinations of this specific KP with each TF and each gene. 
    # KP pert values are repeated automatically #TF times
    DT = data.table(KP_pert[KP == kp], TF_regulons)
    # invert Regulon bool since case should be tested for being greater than control, and case is 0 while control is 1 apparently
    # it's very clear which way to test, the values reveal nothing significant if chosen in the wrong direction.
    DT = DT %>% group_by(TF) %>% summarise(p = wilcox.test(abs(M) ~ !Regulon, alternative="g")$p.value, 
                                           estimate = wilcox.test(M ~ !Regulon, conf.int=TRUE)$estimate)
    
    DT$KP = kp
    write.table(DT, paste0("KP2TF_parts/", kp, ".tsv"), sep="\t", quote=F, row.names=F)
}
