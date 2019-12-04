setwd("/Users/christian/GoogleDrev/PKTFX/data/pertubation")

holstege_PK_KO = read.table("../processed/holstege_2010/PK_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
luscombe_TF_KO = read.table("../processed/luscombe_2010/TF_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
fiedler_PK_KO = read.table("../processed/fiedler_2009/avg_KO.tsv", header=T, row.names=1, check.names=F)
fiedler_PK_KO = fiedler_PK_KO[,colnames(fiedler_PK_KO) != "Gene"]
zelezniak_PK_KO = read.table("../processed/zelezniak_2018/PK_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
chua_TF_KO = read.table("../processed/chua_2006/TF_KO.tsv", header=T, row.names=1, check.names=F)
chua_TF_OE = read.table("../processed/chua_2006/TF_OE.tsv", header=T, row.names=1, check.names=F)
goncalves_KP_KO_TF = read.table("../processed/goncalves_2017/KP_KO_TF.tsv", header=T, row.names=1, sep='\t', check.names=F)
goncalves_KP_KO_KP = read.table("../processed/goncalves_2017/KP_KO_KP.tsv", header=T, row.names=1, sep='\t', check.names=F)

KP_names = c(colnames(holstege_PK_KO), colnames(fiedler_PK_KO), colnames(zelezniak_PK_KO), colnames(goncalves_KP_KO_KP), colnames(goncalves_KP_KO_TF))
TF_names = c(colnames(luscombe_TF_KO), colnames(chua_TF_KO), colnames(chua_TF_OE))
KP_names = sort(unique(unlist(strsplit(KP_names, "_"))))
TF_names = sort(unique(unlist(strsplit(TF_names, "_"))))
# Luscombe have KOs of 5 KPs so we remove those from TF
TF_names = TF_names[!(TF_names %in% KP_names)]
write.table(KP_names, file="KP.txt", row.names=F, col.names=F, quote=F)
write.table(TF_names, file="TF.txt", row.names=F, col.names=F, quote=F)

# merge all pertubation data
dfs = list(holstege_PK_KO, luscombe_TF_KO, fiedler_PK_KO, zelezniak_PK_KO, chua_TF_KO, chua_TF_OE, goncalves_KP_KO_TF, goncalves_KP_KO_KP)
pertubation = dfs[[1]]
for (i in 2:length(dfs)) {
    # by=0 means row names
    pertubation = merge(pertubation, dfs[[i]], by=0, all=T, suffixes=c("", paste0(" ", i)))
    # the merge moves the row names to the first column. we move it back
    row.names(pertubation) = pertubation$Row.names
    pertubation = pertubation[,2:ncol(pertubation)]
}
dim(pertubation)
# duplicate colnames have had .1, .2, etc. appended 
colnames(pertubation) = gsub("\\.[0-9]", "", colnames(pertubation))
colnames(pertubation) = gsub(" [0-9]", "", colnames(pertubation))
dim(pertubation)

X_names = sort(unique(rownames(pertubation)[!(rownames(pertubation) %in% KP_names | rownames(pertubation) %in% TF_names)]))
write.table(X_names, file="X.txt", row.names=F, col.names=F, quote=F)

# notfound
colnames(pertubation)[!(colnames(pertubation) %in% rownames(pertubation))]
# they are either duplicates or multi KO so all is good

PTX_names = c(KP_names, TF_names, X_names)
write.table(PTX_names, file="PTX.txt", row.names=F, col.names=F, quote=F)

# sort rows according to PTX
pertubation = pertubation[order(match(rownames(pertubation), PTX_names)),]
pertubation = pertubation[,order(match(colnames(pertubation), PTX_names))]
# colnames made unique again
colnames(pertubation) = gsub("\\.[0-9]", "", colnames(pertubation))

write.table(pertubation, "logFC.mat", sep=" ", quote=F)


# make idx of what is KOed/OEed
library(Matrix)
i = c(); j = c()
for (column in 1:ncol(pertubation)) {
    KOOE_names = strsplit(colnames(pertubation)[column], "_")[[1]]
    new_i = match(KOOE_names, rownames(pertubation))
    i = c(i, new_i)
    j = c(j, rep(column, length(new_i)))
}
pertubation_indices = sparseMatrix(i, j, x=1, dims=dim(pertubation), dimnames=dimnames(pertubation), use.last.ij=T)
write.table(as.matrix(pertubation_indices), file="pertubation_indices.mat", sep=" ", quote=F)


# make a mask for WT edges
T_edges = read.table("../goldstandard/T/edges.tsv", header=T, sep="\t", quote="")
i = c(); j = c()
for (T_row in 1:nrow(T_edges)) {
    i = c(i, match(T_edges$Target[T_row], PTX_names))
    j = c(j, match(T_edges$TF[T_row], TF_names))
}
# there's not a perfect overlap between the sets so NAs are made
notfound = is.na(i) | is.na(j)
i = i[!notfound]
j = j[!notfound]
WT_prior = sparseMatrix(i, j, x=1, dims=c(length(PTX_names), length(TF_names)), dimnames=list(PTX_names, TF_names), use.last.ij=T)
write.table(as.matrix(WT_prior), file="WT_prior.mat", sep=" ", quote=F)

# this many prior edges
sum(WT_prior)
# how many per TF?
edge_per_TF = apply(WT_prior, 2, sum)
plot(density(edge_per_TF))
mean(edge_per_TF)





