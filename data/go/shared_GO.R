
# count number of GO terms in common for each pair of genes

library(data.table)
library(ggplot2)

flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
melt_ = function(x) melt(x, id.vars="ORF", variable.name="KP", value.name="shared")

KPs = read.vector("~/cwd/data/network/KP.txt")
TFs = read.vector("~/cwd/data/network/TF.txt")
Vs = read.vector("~/cwd/data/network/V.txt")
GO = fread("~/cwd/data/processed/SGD/GO.tsv", 
           sep="\t", col.names=c("ORF", "Aspect", "GOID"), key="ORF")
GO$GOID = as.numeric(gsub("GO:", "", GO$GOID))
GO = GO[Vs,]
GO[,ORF:=factor(GO$ORF, levels=Vs)]
Vs = factor(Vs, levels=Vs)
nV = length(Vs)

GO_C = GO[Aspect=="C",]
GO_F = GO[Aspect=="F",]
GO_P = GO[Aspect=="P",]

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

# takes time, load from compressed result files
# shared_GO_C = get_shared_GO(GO_C)
# shared_GO_F = get_shared_GO(GO_F)
# shared_GO_P = get_shared_GO(GO_P)
# fwrite(shared_GO_C, "shared_GO_C.csv", sep=",")
# fwrite(shared_GO_F, "shared_GO_F.csv", sep=",")
# fwrite(shared_GO_P, "shared_GO_P.csv", sep=",")
shared_GO_C = fread("shared_GO_C.csv.gz", sep=",", key="ORF")
shared_GO_F = fread("shared_GO_F.csv.gz", sep=",", key="ORF")
shared_GO_P = fread("shared_GO_P.csv.gz", sep=",", key="ORF")

shared_GO = rbind(
    melt_(shared_GO_C[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("C","KP")],
    melt_(shared_GO_F[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("F","KP")],
    melt_(shared_GO_P[KPs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("P","KP")],
    melt_(shared_GO_C[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("C","TF")],
    melt_(shared_GO_F[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("F","TF")],
    melt_(shared_GO_P[TFs,c("ORF",..KPs)])[,c("Aspect","Substrate"):=.("P","TF")])

# no self loops
shared_GO = shared_GO[ORF!=KP,]

fwrite(shared_GO, "shared_GO.tsv", sep="\t")




