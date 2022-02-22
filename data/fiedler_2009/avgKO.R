#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data/fiedler_2009"))
# take average of KO values since there are not enough values to call it replicates.
KO_table = read.table("ctk1_cak1_KO.tsv", sep="\t", header=T, quote="")
avg_table = aggregate(. ~ ORF + Gene, data=KO_table, mean)
write.table(avg_table[,2:3], file="avg_KO.tsv", sep="\t", quote=FALSE, row.names=FALSE)

