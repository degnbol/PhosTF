
library(tidyverse)
library(readxl)

setwd("~/cwd/data/raw/ptacek_2005")
filename = "41586_2005_BFnature04187_MOESM4_ESM.xls"

sheet_names = readxl::excel_sheets(filename)

edges = tibble()
for (sheet_name in sheet_names) {
    sheet = readxl::read_excel(filename, sheet=sheet_name, range=cell_cols("A"))
    colnames(sheet) = "Target"
    sheet$KP = sheet_name
    edges = rbind(edges, sheet)
}

edges = edges[,2:1]

edges$KP = gsub("^CDC28", "CDC28 ", edges$KP)
edges$KP = gsub("^PHO85", "PHO85 ", edges$KP)
edges$KP = gsub(" alone", "", edges$KP)

write.table(edges, "~/cwd/data/processed/ptacek_2005/edges.tsv", sep="\t", quote=F, row.names=F)






