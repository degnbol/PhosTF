
library(reshape2)

setwd("~/cwd/data/processed/zelezniak_2018")

dataset = read.table("proteins_dataset.data_prep.tsv", header=T, sep="\t")
WT = aggregate(value ~ KO_ORF + ORF, data=dataset[dataset$KO_ORF == "WT",], mean)
measurements = dataset[dataset$KO_ORF != "WT", !(colnames(dataset) %in% c("sample", "KO_gene_name"))]
measurements2D = acast(measurements, ORF ~ KO_ORF + replicate)
dim(measurements2D)
# allow replicates to have identical column name
colnames(measurements2D) = gsub("_.*", "", colnames(measurements2D))

# sanity check
all(rownames(measurements2D) == WT$ORF)

logFC = log2(measurements2D / WT$value)
write.table(logFC, file="PK_KO.tsv", sep="\t", quote=F)


dataset$ORF == dataset$KO_ORF

