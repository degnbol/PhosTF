#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# for testing:
file = "../../results/4-yeastNetworkReconstruction/WT_FDR20_sign_mask.ssv"

# defaults
rowName = "row"
colName = "column"
valueName = "value"

# read args
args = commandArgs(TRUE)
file = args[1]
if(length(args) > 1) {
    rowName = args[2]
    if(length(args) > 2) {
        colName = args[3]
        if(length(args) > 3) {
            valueName = args[4]
        }
    }
}
        
        
# if first column name is empty to indicate row names then fread detects it correctly but complains.
dt = suppressWarnings(fread(file))

stopifnot(names(dt)[1] == "V1")
setnames(dt, "V1", rowName)

dt.melt = melt(dt, rowName, variable.name=colName, value.name=valueName)

fwrite(dt.melt, sep='\t')

