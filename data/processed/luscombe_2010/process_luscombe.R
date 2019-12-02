setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/luscombe_2010")

TFKO = t(read.table("TF_KO_transpose.tsv", row.names=1, header=T))
TFKO_pval = t(read.table("TF_KO_pval_transpose.tsv", row.names=1, header=T))
# sanity check
dim(TFKO)
dim(TFKO_pval)

# convert name using uniprot
gene2ORF_table = read.table("../uniprot/gene2ORF.tsv", header=F, sep="\t", quote="")
gene2ORF = function(x) gene2ORF_table[match(x,gene2ORF_table[,1]),2]

# first a sanity check
all(colnames(TFKO) == colnames(TFKO_pval))

names = toupper(colnames(TFKO))
ORFs = gene2ORF(names)
ORFs[is.na(ORFs)] = names[is.na(ORFs)]
# update colnames
colnames(TFKO) = ORFs
colnames(TFKO_pval) = ORFs

write.table(TFKO, file="TF_KO.tsv", quote=F, sep="\t", row.names=T)
write.table(TFKO_pval, file="TF_KO_pval.tsv", quote=F, sep="\t", row.names=T)
