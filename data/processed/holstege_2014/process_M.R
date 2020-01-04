
library(data.table)

flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/processed/holstege_2014")

M_table = fread("M.tsv", sep="\t", header=T, quote="")

# some genes are measured multiple times
dups_index = duplicated(M_table$ORF)
dups = M_table$ORF[dups_index]
length(dups) == length(unique(dups))
# only one of each, easy to mean
for (dup in dups) {
    dup_values = M_table[M_table$ORF == dup, 2:ncol(M_table)]
    M_table[M_table$ORF == dup, 2:ncol(M_table)] = (dup_values[1,] + dup_values[2,]) / 2
}
# remove dups
M_table = M_table[!dups_index,]


write.table(M_table, "KO.tsv", sep="\t", quote=F, row.names=F)
