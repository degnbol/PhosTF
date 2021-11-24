#!/usr/bin/env Rscript
library(tidyverse)
library(readxl)

filename = "raw/41586_2005_BFnature04187_MOESM4_ESM.xls"

sheet_names = readxl::excel_sheets(filename)

edges = tibble()
for (sheet_name in sheet_names) {
    sheet = readxl::read_excel(filename, sheet=sheet_name, range=cell_cols("A"))
    colnames(sheet) = "Target"
    sheet$KP = sheet_name
    edges = rbind(edges, sheet)
}

edges$KP2 = ""
for (paired in c("CDC28", "PHO85")) {
    regex = paste0("^",paired)
    idx = grepl(regex, edges$KP)
    edges$KP2[idx] = gsub(regex, "", edges$KP[idx])
    edges$KP[idx] = paired
}
edges$KP2[edges$KP2 == "alone"] = ""

edges = edges[,c("KP", "KP2", "Target")]

write.table(edges, "edges_pop.tsv", sep="\t", quote=F, row.names=F)






