setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/holstege_2010")
# read
PKKO = read.table("PK_KO_Mpval.tsv", header=T, sep="\t", quote="")
# only logFC values and ORF
PKKOM = PKKO[,seq(1,ncol(PKKO),2)]
# what are the KOed genes?
KO_genes = toupper(gsub("\\..*", "", names(PKKOM)))[2:ncol(PKKOM)]
# do we have an ORF for each?
KO_genes[!(KO_genes %in% PKKO$Gene)]
# we are missing ORF for "CDK8" and "SHA3", look it up manually on yeastgenome.org
# https://www.yeastgenome.org/locus/S000005963 CDK8 = YPL042C
# https://www.yeastgenome.org/locus/S000005947 SHA3 = YPL026C
KO_genes[KO_genes == "CDK8"] = "YPL042C"
KO_genes[KO_genes == "SHA3"] = "YPL026C"
# copy ORFs from the translations among measured genes
KO_ORFs = PKKO$ORF[match(KO_genes, PKKO$Gene)]
# those that weren't found are in the initial list
KO_ORFs[is.na(KO_ORFs)] = KO_genes[is.na(KO_ORFs)]
# check that there are no na
any(is.na(KO_ORFs))
# update KO names to ORF
names(PKKOM)[2:ncol(PKKOM)] = as.character(KO_ORFs)
write.table(PKKOM, "PK_KO.tsv", quote=F, sep="\t", row.names=F)
