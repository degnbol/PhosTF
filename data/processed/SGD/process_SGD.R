
setwd("~/cwd/data/processed/SGD")

ORF_gene_aliases = read.table("ORF_gene_aliases.tsv", sep="\t", check.names=F, quote="", col.names=c("ORF", "Gene", "Aliases"), stringsAsFactors=F)

aliases = strsplit(ORF_gene_aliases$Aliases, '|', fixed=T)

convert.from = c()
convert.to   = c()

for (i in 1:nrow(ORF_gene_aliases))
    for (alias in aliases[[i]]) {
        convert.from = c(convert.from, alias)
        convert.to   = c(convert.to, ORF_gene_aliases$ORF[[i]])
    }

conversions = data.frame(ORF=convert.to, Alias=convert.from)
write.table(conversions, file="alias2ORF.tsv", sep="\t", quote=F, row.names=F)
