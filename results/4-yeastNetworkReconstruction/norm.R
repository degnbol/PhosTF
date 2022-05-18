#!/usr/bin/env Rscript
# Normalize matrix in each column.
# USE: ./norm.R INFILE.csv | sed '/1/s/""//' | gzip > OUTFILE.csv.gz
# sed to replace empty rownames header '""' with ''.
suppressPackageStartupMessages(library(data.table))

debug = FALSE
if(debug) {
    dt = fread("logFC_inner-strict.csv")
} else {
    dt = fread(commandArgs(TRUE))
}

# assume format first column is rownames without column name (becomes V1)
dt.melt = melt(dt)
dt.melt[, norm:=(value - mean(value)) / sd(value), by=variable]
dt.cast = dcast(dt.melt, V1 ~ variable, value.var="norm")
setnames(dt.cast, "V1", "")

fwrite(dt.cast)

