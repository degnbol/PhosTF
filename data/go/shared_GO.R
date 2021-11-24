#!/usr/bin/env Rscript
# count number of GO terms in common for each pair of genes

SuppressPackageStartupMessages(library(data.table))
library(ggplot2)
SuppressPackageStartupMessages(library(here))

cwd = function(s) paste0(here(), "/", s)
flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(cwd(x)))
melt_ = function(x) melt(x, id.vars="ORF", variable.name="KP", value.name="shared")
prepare = function(GO_table) {
    GO_table$GOID = as.numeric(gsub("GO:", "", GO_table$GOID))
    GO_table = GO_table[Vs,]
    GO_table[,ORF:=factor(GO_table$ORF, levels=Vs)]
    GO_table
}
get_shared_GO = function(GO_table) {
    pairs = data.table(ORF=Vs, key="ORF")
    for(v in Vs) {
        # GO IDs for an ORF
        GOIDs = GO_table[ORF==v,GOID]
        # count how many of those GO IDs are found for all other ORFs
        shared = GO_table[GOID%in%GOIDs, .N, by=ORF]
        v = as.character(v)
        # left join on ORF,
        # which means NA are inserted for each pair where no GO IDs are in common
        pairs[shared, on="ORF", (v):=N]
    }
    pairs[is.na(pairs)] = 0
    pairs
}


KPs = read.vector("data/network/KP.txt")
TFs = read.vector("data/network/TF.txt")
Vs = read.vector("data/network/V.txt")
GO = fread(cwd("data/processed/SGD/GO.tsv"), 
           sep="\t", col.names=c("ORF", "Aspect", "GOID"), key="ORF")
GOP250 = fread(cwd("data/processed/SGD/GOP250.tsv"), sep="\t", key="ORF")
GOP400 = fread(cwd("data/processed/SGD/GOP400.tsv"), sep="\t", key="ORF")
GOP500 = fread(cwd("data/processed/SGD/GOP500.tsv"), sep="\t", key="ORF")

GO = prepare(GO)
GOP250 = prepare(GOP250)
GOP400 = prepare(GOP400)
GOP500 = prepare(GOP500)

Vs = factor(Vs, levels=Vs)
nV = length(Vs)


GO_C = GO[Aspect=="C",]
GO_F = GO[Aspect=="F",]
GO_P = GO[Aspect=="P",]

# takes time, load from compressed result files
# shared_GO_C = get_shared_GO(GO_C)
# shared_GO_F = get_shared_GO(GO_F)
# shared_GO_P = get_shared_GO(GO_P)
# shared_GOP250 = get_shared_GO(GOP250)
# shared_GOP400 = get_shared_GO(GOP400)
# shared_GOP500 = get_shared_GO(GOP500)
# fwrite(shared_GO_C, "shared_GO_C.csv", sep=",")
# fwrite(shared_GO_F, "shared_GO_F.csv", sep=",")
# fwrite(shared_GO_P, "shared_GO_P.csv", sep=",")
# fwrite(shared_GOP250, "shared_GOP250.csv", sep=",")
# fwrite(shared_GOP400, "shared_GOP400.csv", sep=",")
# fwrite(shared_GOP500, "shared_GOP500.csv", sep=",")
shared_GO_C = fread("shared_GO_C.csv.gz", sep=",", key="ORF")
shared_GO_F = fread("shared_GO_F.csv.gz", sep=",", key="ORF")
shared_GO_P = fread("shared_GO_P.csv.gz", sep=",", key="ORF")
shared_GOP250 = fread("shared_GOP250.csv.gz", sep=",", key="ORF")
shared_GOP400 = fread("shared_GOP400.csv.gz", sep=",", key="ORF")
shared_GOP500 = fread("shared_GOP500.csv.gz", sep=",", key="ORF")

shared_GO = rbind(
    melt_(shared_GO_C[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("C","KP")],
    melt_(shared_GO_F[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("F","KP")],
    melt_(shared_GO_P[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("P","KP")],
    melt_(shared_GO_C[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("C","TF")],
    melt_(shared_GO_F[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("F","TF")],
    melt_(shared_GO_P[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("P","TF")])


GO_melt = function(GO_table) {
    rbind(
        melt_(GO_table[KPs,c("ORF",..KPs)])[,Substrate:="KP"],
        melt_(GO_table[TFs,c("ORF",..KPs)])[,Substrate:="TF"])
}

shared_GOP250_melt = GO_melt(shared_GOP250)
shared_GOP400_melt = GO_melt(shared_GOP400)
shared_GOP500_melt = GO_melt(shared_GOP500)

# no self loops
shared_GO = shared_GO[ORF!=KP,]
shared_GOP250_melt = shared_GOP250_melt[ORF!=KP,]
shared_GOP400_melt = shared_GOP400_melt[ORF!=KP,]
shared_GOP500_melt = shared_GOP500_melt[ORF!=KP,]

fwrite(shared_GO, "shared_GO.tsv", sep="\t")
fwrite(shared_GOP250_melt, "shared_GOP250.tsv", sep="\t")
fwrite(shared_GOP400_melt, "shared_GOP400.tsv", sep="\t")
fwrite(shared_GOP500_melt, "shared_GOP500.tsv", sep="\t")




