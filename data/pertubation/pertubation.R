
## packages
library(Matrix)
library(reshape2)

## functions

# make logical matrix with 1 indicating that a row name is in a column name.
# column names are split with _ so we assume column names with multiple names in them are written that way
column_row_map = function(x) {
    i = c(); j = c()
    for (column in 1:ncol(x)) {
        column_names = strsplit(colnames(x)[column], "_")[[1]]
        new_i = match(column_names, rownames(x))
        i = c(i, new_i)
        j = c(j, rep(column, length(new_i)))
    }
    sparseMatrix(i, j, x=1, dims=dim(x), dimnames=dimnames(x), use.last.ij=T)
}
# also split with space so only the first word is used
column_row_map_space = function(x) {
    i = c(); j = c()
    for (column in 1:ncol(x)) {
        column_names = strsplit(gsub(" .*", "", colnames(x)[column]), "_")[[1]]
        new_i = match(column_names, rownames(x))
        i = c(i, new_i)
        j = c(j, rep(column, length(new_i)))
    }
    sparseMatrix(i, j, x=1, dims=dim(x), dimnames=dimnames(x), use.last.ij=T)
}

# get values where column name == row name
# splitting column name by _
column_row_values = function(x) {
    values = c()
    for (j in 1:ncol(x)) {
        column_names = strsplit(colnames(x)[j], "_")[[1]]
        for (i in match(column_names, rownames(x))) {
            values = c(values, x[i,j])
        }
    }
    values
}

# make a logical matrix with the given row names and column names where 1 indicates a match from the given row_matches and column_matches
match_columnrow = function(row_names, column_names, row_matches, column_matches) {
    i = c(); j = c()
    for (row_match in row_matches) {
        i = c(i, match(row_match, row_names))
    }
    for (column_match in column_matches) {
        j = c(j, match(column_match, column_names))
    }
    # if there's not a perfect overlap between the sets then NAs are made
    notfound = is.na(i) | is.na(j)
    i = i[!notfound]
    j = j[!notfound]
    sparseMatrix(i, j, x=1, dims=c(length(row_names), length(column_names)), dimnames=list(row_names, column_names), use.last.ij=T)
}

# sort rows and columns in a table given a collection of elements that are matched against which are in the correct order
sort_table = function(x, sort_order) {
    # sort rows according to a collection of names that are in the intended order
    x = x[order(match(rownames(x), sort_order)),]
    x = x[,order(match(colnames(x), sort_order))]
    # colnames made unique again
    colnames(x) = gsub("\\.[0-9]", "", colnames(x))
    x
}

melt_table = function(x) {
    out = melt(as.matrix(x))
    colnames(out) = c("rownames", "colnames", "value")
    out
}

melt_tables = function(dfs) {
    out = melt_table(dfs[[1]])
    for (i in 2:length(dfs)) {
        out = rbind(out, melt_table(dfs[[i]]))
    }
    out
}

# dfs is a list of data.frames
merge_full = function(dfs) {
    out = dfs[[1]]
    for (i in 2:length(dfs)) {
        # by=0 means row names
        out = merge(out, dfs[[i]], by=0, all=T, suffixes=c("", paste0(" ", i)))
        # the merge moves the row names to the first column. we move it back
        row.names(out) = out$Row.names
        out = out[,2:ncol(out)]
    }
    # duplicate colnames have had .1, .2, etc. appended 
    colnames(out) = gsub("\\.[0-9]", "", colnames(out))
    colnames(out) = gsub(" [0-9]", "", colnames(out))
    out
}

# replace NA in table x with average of the values in other columns with identical column names from the same row
NA2avg = function(x) {
    # mean wants matrix
    out = as.matrix(x)
    for(i in 1:nrow(x))
        for(j in 1:ncol(x))
            if(is.na(out[i,j])) {
                out[i,j] = mean(x[i, colnames(x) == colnames(x)[j]], na.rm=T)
            }
    out
}
# same purpose as above but we have another table with averages to copy from
# assume rownames match
NA2avg = function(x, avg_table) {
    stopifnot(all(rownames(x) == rownames(avg_table)))
    # copy columns from avg_table into a new table where rownames and colnames matches x
    avg_x = avg_table[,colnames(x)]
    out = x
    out[is.na(x)] = avg_x[is.na(x)]
    out
}


setwd("/Users/christian/GoogleDrev/PKTFX/data/pertubation")

## read
holstege_PK_KO = read.table("../processed/holstege_2010/PK_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
luscombe_TF_KO = read.table("../processed/luscombe_2010/TF_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
fiedler_PK_KO = read.table("../processed/fiedler_2009/avg_KO.tsv", header=T, row.names=1, check.names=F)
fiedler_PK_KO = fiedler_PK_KO[,colnames(fiedler_PK_KO) != "Gene"]
zelezniak_PK_KO = read.table("../processed/zelezniak_2018/PK_KO.tsv", header=T, row.names=1, sep='\t', check.names=F)
chua_TF_KO = read.table("../processed/chua_2006/TF_KO.tsv", header=T, row.names=1, check.names=F)
chua_TF_OE = read.table("../processed/chua_2006/TF_OE.tsv", header=T, row.names=1, check.names=F)
goncalves_KP_KO_TF = read.table("../processed/goncalves_2017/KP_KO_TF.tsv", header=T, row.names=1, sep='\t', check.names=F)
goncalves_KP_KO_KP = read.table("../processed/goncalves_2017/KP_KO_KP.tsv", header=T, row.names=1, sep='\t', check.names=F)
data_frames = list(holstege_PK_KO, luscombe_TF_KO, fiedler_PK_KO, zelezniak_PK_KO, chua_TF_KO, chua_TF_OE, goncalves_KP_KO_TF, goncalves_KP_KO_KP)

KP_names = c(colnames(holstege_PK_KO), colnames(fiedler_PK_KO), colnames(zelezniak_PK_KO), colnames(goncalves_KP_KO_KP), colnames(goncalves_KP_KO_TF))
TF_names = c(colnames(luscombe_TF_KO), colnames(chua_TF_KO), colnames(chua_TF_OE))
KP_names = sort(unique(unlist(strsplit(KP_names, "_"))))
TF_names = sort(unique(unlist(strsplit(TF_names, "_"))))
# Luscombe have KOs of 5 KPs so we remove those from TF
TF_names = TF_names[!(TF_names %in% KP_names)]
write.table(KP_names, file="KP.txt", row.names=F, col.names=F, quote=F)
write.table(TF_names, file="TF.txt", row.names=F, col.names=F, quote=F)

colnames(chua_TF_OE) = paste(colnames(chua_TF_OE), "OE")  # make sure to separate OE from KO

# all values in long format
data_frames_melt = melt_tables(data_frames)
data_frames_melt_avg = aggregate(value ~ rownames + colnames, data=data_frames_melt, mean)
data_frames_cast = acast(data_frames_melt_avg, rownames ~ colnames)

# merge all pertubation data
pertubation = merge_full(data_frames)

X_names = sort(unique(rownames(pertubation)))
X_names = X_names[!((X_names %in% KP_names) | (X_names %in% TF_names))]
write.table(X_names, file="X.txt", row.names=F, col.names=F, quote=F)

# check if any pertubation gene is not measured
colnames(pertubation)[!(gsub(" OE", "", colnames(pertubation)) %in% rownames(pertubation))]
# they are only multi KO so all is good

PTX_names = c(KP_names, TF_names, X_names)
write.table(PTX_names, file="PTX.txt", row.names=F, col.names=F, quote=F)

# sort rows according to PTX
pertubation = sort_table(pertubation, PTX_names)
pertubation_inner = sort_table(data_frames_cast, PTX_names)

# make idx of what is KOed/OEed
KO_indices = as.matrix(column_row_map(pertubation))
KOOE_indices = as.matrix(column_row_map_space(pertubation))
KO_inner_indices = as.matrix(column_row_map(pertubation_inner))
KOOE_inner_indices = as.matrix(column_row_map_space(pertubation_inner))
write.table(KOOE_indices, file="KOOE_indices.ssv", sep=" ", quote=F)
write.table(KOOE_inner_indices, file="KOOE_inner_indices.ssv", sep=" ", quote=F)
# we add NaN entries to J, which will result in their values being ignored in the SSE calculation during gradient descent
J = KOOE_indices
J_inner = KOOE_inner_indices
J[is.na(pertubation)] = 1
J_inner[is.na(pertubation_inner)] = 1
write.table(J, file="J.ssv", sep=" ", quote=F)
write.table(J_inner, file="J_inner.ssv", sep=" ", quote=F)

# make them logical
KO_indices = KO_indices == 1
KOOE_indices = KOOE_indices == 1
KO_inner_indices = KO_inner_indices == 1
KOOE_inner_indices = KOOE_inner_indices == 1

sum(is.na(pertubation))
# we have to get rid of NaNs. We can try to copy logFC for a node under the same exp conditions
pertubation = NA2avg(pertubation, pertubation_inner)
# the remaining NA values are set to zero.
sum(is.na(pertubation))
sum(is.na(pertubation_inner))
# correct remaining nodes to zero
pertubation[is.na(pertubation)] = 0
pertubation_inner[is.na(pertubation_inner)] = 0
# reducing KOs with 6
pertubation[KO_indices] = pertubation[KO_indices] - 6
pertubation_inner[KO_inner_indices] = pertubation_inner[KO_inner_indices] - 6
# write
write.table(pertubation, file="logFC.ssv", sep=" ", quote=F)
write.table(pertubation_inner, "logFC_inner.ssv", sep=" ", quote=F)


# make a mask for WT edges
T_edges = read.table("../goldstandard/T/edges.tsv", header=T, sep="\t", quote="")
WT_prior = match_columnrow(PTX_names, TF_names, T_edges$Target, T_edges$TF)
write.table(as.matrix(WT_prior), file="WT_prior.ssv", sep=" ", quote=F)
# this many prior edges
sum(WT_prior)
# how many per TF?
edge_per_TF = apply(WT_prior, 2, sum)
plot(density(edge_per_TF))
mean(edge_per_TF)






