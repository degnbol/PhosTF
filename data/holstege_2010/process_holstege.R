
setwd("~/cwd/data/processed/holstege_2010")
# read
PKKO = read.table("PK_KO_Mpval.tsv", header=T, sep="\t", quote="", stringsAsFactors=F)
# only logFC values and ORF
stopifnot(all(grep("\\.M", colnames(PKKO)) == seq(3,ncol(PKKO),2)))  # make sure M values are actually in every other row
PKKOM = PKKO[,seq(1,ncol(PKKO),2)]
# what are the KOed genes?
KO_genes = toupper(gsub("\\.M", "", names(PKKOM)))[2:ncol(PKKOM)]
# some are multiple KOs at the same time
KO_genes_mult = strsplit(KO_genes, "\\.")

# do we have an ORF for each?
KO_genes[!(KO_genes %in% PKKO$Gene)]
# we are missing ORF for "CDK8" and "SHA3", look it up manually on yeastgenome.org
# https://www.yeastgenome.org/locus/S000005963 CDK8 = YPL042C
# https://www.yeastgenome.org/locus/S000005947 SHA3 = YPL026C
KO_genes[KO_genes == "CDK8"] = "YPL042C"
KO_genes[KO_genes == "SHA3"] = "YPL026C"
# copy ORFs from the translations among measured genes
KO_ORFs = c()
for (gs in KO_genes_mult) {
    KO_ORFs = c(KO_ORFs, paste(PKKO$ORF[match(gs, PKKO$Gene)], collapse="_"))
}
# those that weren't found are in the initial list
KO_ORFs[KO_ORFs == "NA"] = KO_genes[KO_ORFs == "NA"]
# check that there are no na
any(is.na(KO_ORFs))
# update KO names to ORF
names(PKKOM)[2:ncol(PKKOM)] = KO_ORFs
# average value of the few replicates in this data
agg = aggregate(. ~ ORF, data=PKKOM, mean)
write.table(agg, "PK_KO.tsv", quote=F, sep="\t", row.names=F)
