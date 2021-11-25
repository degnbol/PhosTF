#!/usr/bin/env Rscript
# Alternative download method, although it seems only unsigned.

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("DREAM4")
library(DREAM4)
suppressPackageStartupMessages(library(data.table))

for(nNodes in c(10, 100)) { 
    for(i in 1:5) {
        name = sprintf("dream4_%03d_%02d", nNodes, i)
        data(list=name)
        adj = metadata(get(name))[[1]]
        fwrite(adj, paste0(name, ".adj"), sep=" ", col.names=FALSE)
    }
}

