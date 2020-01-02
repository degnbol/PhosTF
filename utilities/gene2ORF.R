#!/usr/bin/env Rscript

# functions
simplify = function(x) gsub(" ", "", gsub("_", "", toupper(trimws(x))))
# same as match but ignore case and _ and " "
match_ignore = function(x, table) match(simplify(x), simplify(table))
# assuming ORF in column 1, and Gene in column 2
convert_column = function(column, ignore=F) {
    if (ignore) return(converter[match_ignore(column, converter[,2]),1])
    else return(converter[match(column, converter[,2]),1])
}
convert_table = function(gene.table, ignore=F) {
    ORF.table = gene.table
    for (i in 1:ncol(ORF.table)) {
        ORF.table[,i] = convert_column(gene.table[,i], ignore=ignore)
        # we use the original value when no conversion found
        no_conversion_idx = is.na(ORF.table[,i])
        ORF.table[no_conversion_idx,i] = gene.table[no_conversion_idx,i]
    }
    ORF.table
}

# commands
args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3)
# doesn't matter if file has header or not
infname = args[1]
# doesn't matter if file has header or not
outfname = args[2]
# doesn't matter if file has header or not
# should have ORF in column 1 and Gene in column 2
converter_fname = args[3]

# main
gene.table = read.table(infname, sep="\t", quote="", stringsAsFactors=F, colClasses="character")
converter = read.table(converter_fname, sep="\t", quote="", stringsAsFactors=F)
ORF.table = convert_table(convert_table(gene.table, ignore=F), ignore=T)
write.table(ORF.table, file=outfname, row.names=F, col.names=F, quote=F, sep="\t")
