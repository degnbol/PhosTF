#!/usr/bin/env Rscript
# list all evidence together

library(data.table)

# functions
unwhich = function(which, dim=max(which)) {
    y = array(logical(dim), dim=dim)
    y[which] = TRUE
    y
}

read_go = function(path, label) {
    go = fread(path, col.names=c("Evidence", "ORF"))
    go[,Label:=label]
    return(go)
}

setwd("~/cwd/data/processed/amigo2")

# read GO

GO = rbind(
    read_go("transcription_factors.tsv", "TF"),
    read_go("activators.tsv", "TFA"),
    read_go("repressors.tsv", "TFR"),
    read_go("protein_kinase.tsv", "PK"),
    read_go("protein_phosphatase.tsv", "PP"),
    read_go("protein_phosphorylation.tsv", "PK"),
    read_go("protein_dephosphorylation.tsv", "PP")
)

# split into 3 categories of evidence
# http://wiki.geneontology.org/index.php/Guide_to_GO_Evidence_Codes
experimental_evidence = c("HMP", "IMP", "HDA", "IGI", "IDA", "IPI")
computational_evidence = c("IEA", "IGC", "RCA", "IBA", "ISO", "ISA", "ISS", "ISM")
statement_evidence = c("IC", "NAS", "TAS")

GO[Evidence%in%experimental_evidence,Experimental:=Evidence]
GO[,Evidence:=NULL]
for(label in unique(GO$Label)) {
    GO[Label==label,paste(label,"amigo",sep="_"):=Experimental]
}
GO[,Experimental:=NULL]
# GO[,Nonexperimental:=NULL]
GO[,Label:=NULL]
GO[is.na(GO)] = ""

GO = GO[, lapply(.SD, function(x) paste(x[x!=""], collapse=",")), by=ORF]


fwrite(GO, "amigo_evidences.tsv", sep="\t", quote=F)







