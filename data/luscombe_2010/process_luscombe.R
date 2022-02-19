#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(here))
setwd(paste0(here(), "/data/luscombe_2010"))

# IMPORTANT that R does not check names since it will rename - to .
KO = t(read.table("KO_transpose.tsv", row.names=1, header=T, check.names=F, stringsAsFactors=F))
KO_pval = t(read.table("KO_pval_transpose.tsv", row.names=1, header=T, check.names=F, stringsAsFactors=F))
# sanity checks
dim(KO); dim(KO_pval)
all(colnames(KO) == colnames(KO_pval))

# convert name
gene2ORF_table = read.table("../../name_conversion/gene2ORF.tsv", header=F, sep="\t", quote="", stringsAsFactors=F)
gene2ORF = function(x) gene2ORF_table[match(x,gene2ORF_table[,2]),1]

names = toupper(colnames(KO))
ORFs = gene2ORF(names)
ORFs[is.na(ORFs)] = names[is.na(ORFs)]
# update colnames
colnames(KO) = ORFs
colnames(KO_pval) = ORFs

# there some MAPK kinases included
kinase_cols = ORFs %in% c("YDR477W", "YLR113W", "YPL042C", "YGR040W", "YPR054W")

write.table(KO[,kinase_cols], "PK_KO.tsv", quote=F, sep="\t", row.names=T)
write.table(KO_pval[,kinase_cols], "PK_KO_pval.tsv", quote=F, sep="\t", row.names=T)
write.table(KO[,!kinase_cols], "TF_KO.tsv", quote=F, sep="\t", row.names=T)
write.table(KO_pval[,!kinase_cols], "TF_KO_pval.tsv", quote=F, sep="\t", row.names=T)

# now go and add ORF to the beginning of file to match standard design of files
