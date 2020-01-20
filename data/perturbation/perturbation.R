#!/usr/bin/env Rscript
## packages
library(data.table)
library(dplyr)
library(reshape2)
library(Matrix)
library(ggplot2)

## functions

flatten = function(x) as.vector(as.matrix(x))

rootname = function(names) {gsub("_.*", "", gsub(" .*", "", names))}

colname_fix = function(x) {
    # duplicate colnames have had .1, .2, etc. appended 
    colnames(x) = gsub("\\.[0-9]+", "", colnames(x))
    colnames(x) = gsub(" [0-9]+", "", colnames(x))
    x
}
# dts is a list of data.tables
merge_outer = function(dts) {
    out = dts[[1]]
    for (i in 2:length(dts)) {
        # make sure col names are unique before a merge
        colnames(out) = make.unique(colnames(out), sep=" ")
        colnames(dts[[i]]) = make.unique(colnames(dts[[i]]), sep=" ")
        out = merge(out, dts[[i]], by="ORF", all=T, suffixes=c("", paste0(" ", i)))
    }
    colname_fix(out)
}


# make logical matrix with 1 indicating that a row name is in a column name.
# column names are split with _ so we assume column names with multiple names in them are written that way
# also split with space so only the first word is used
column_row_map = function(x, column_names) {
    # row_names assumed to be in column 1 SO IT IS SKIPPED here
    row_names = x[,1][[1]]
    x = x[,2:ncol(x)]
    is = c(); js = c()
    for (column in 1:ncol(x)) {
        i = match(column_names[[column]], row_names)
        i = i[!is.na(i)] # NA means there was no match
        is = c(is, i)
        js = c(js, rep(column, length(i)))
    }
    as.matrix(sparseMatrix(is, js, x=1, dims=dim(x), dimnames=list(row_names, colnames(x)), use.last.ij=T))
}
column_row_map_space = function(x) {column_row_map(x, strsplit(gsub(" .*", "", colnames(x)[2:ncol(x)]), "_"))}
column_row_map_nospace = function(x) {column_row_map(x, strsplit(colnames(x)[2:ncol(x)], "_"))}

# copy values from one table to another to replaces NAs. Copies from entries with matching row and column names. assumes matrix type.
NA2avg = function(x, copy_from) {
    # make sure rownames match, assume they are first column
    stopifnot(all(rownames(x) == rownames(copy_from)))
    # copy columns from copy_from into a new table where rownames and colnames matches x
    copy_x = copy_from[,colnames(x)]
    x[is.na(x)] = copy_x[is.na(x)]
    x
}


# main

setwd("~/cwd/data/perturbation")
## read
KP_names = flatten(read.table("../network/KP.txt"))
TF_names = flatten(read.table("../network/TF.txt"))
V_names = flatten(read.table("../network/V.txt"))
O_names = V_names[!(V_names%in%c(KP_names,TF_names))]

fnames = c("../processed/holstege_2010/PK_KO.tsv",
           "../processed/holstege_2014/KO.tsv",
           "../processed/luscombe_2010/TF_KO.tsv",
           "../processed/luscombe_2010/PK_KO.tsv",
           "../processed/fiedler_2009/avg_KO.tsv",
           "../processed/zelezniak_2018/PK_KO.tsv",
           "../processed/goncalves_2017/KP_KO_TF.tsv",
           "../processed/goncalves_2017/KP_KO_KP.tsv",
           "../processed/chua_2006/TF_KO.tsv",
           "../processed/chua_2006/TF_OE.tsv")
pert_tables = list()
for (i in 1:length(fnames)) {pert_tables[[i]] = fread(fnames[i], quote="", header=T, sep='\t')}
OE_fname = grep("OE", fnames)
colnames(pert_tables[[OE_fname]])[2:ncol(pert_tables[[OE_fname]])] = paste(colnames(pert_tables[[OE_fname]])[2:ncol(pert_tables[[OE_fname]])], "OE")  # make sure to separate OE from KO
pert_melt = list()
for (i in 1:length(pert_tables)) {pert_melt[[i]] = melt(pert_tables[[i]], id.vars="ORF", variable.name="Mutant")}
pert_melt = bind_rows(pert_melt)
pert_melt_avg = pert_melt[,list(value=mean(value)),by=list(ORF,Mutant)]
pert_outer = merge_outer(pert_tables) # merge all data, keeping all datapoints
pert_inner = as.data.table(acast(pert_melt_avg, ORF ~ Mutant), "ORF") # merge all data, where there is a single averaged value for each unique entry

# only keep mutation of KP, TF and multi KOs. We want to use _ and OE columns which will be fine, they don't match O_names
pert_outer[,which(colnames(pert_outer) %in% O_names):=NULL]
pert_inner[,which(colnames(pert_inner) %in% O_names):=NULL]

# check if any perturbation gene is not measured
not_measured = colnames(pert_outer)[!(rootname(colnames(pert_outer)) %in% pert_outer$ORF)]
not_measured = not_measured[not_measured != "ORF"] # remove ORF
!any(not_measured %in% c(KP_names, TF_names, O_names))  # they are not in these sets so they are all useless. should be TRUE
pert_outer[,(not_measured):=NULL]  # remove them
pert_inner[,(not_measured):=NULL]  # remove them
# should be TRUE, it is only TFs, which we will assign some arbitrary KO value to let them have influence
all(V_names[!(V_names %in% pert_outer$ORF)] %in% TF_names)

## sorting
# we need to have unique names for sorting
colnames(pert_outer) = make.unique(colnames(pert_outer), sep=" ")
# sort columns according to KP,TF
setcolorder(pert_outer, order(match(rootname(colnames(pert_outer)), V_names), na.last=F))
setcolorder(pert_inner, order(match(rootname(colnames(pert_inner)), V_names), na.last=F))
# sort rows according to KP,TF,O. we use indexing instead of sort to insert 3 new NA lines for TFs that are never measured
setkey(pert_outer,ORF)
setkey(pert_inner,ORF)
pert_outer = pert_outer[V_names,]
pert_inner = pert_inner[V_names,]

# stop the uniqueness of col names
pert_outer = colname_fix(pert_outer)

# make idx of what is KOed/OEed
KO_outer = column_row_map_nospace(pert_outer)
KO_inner = column_row_map_nospace(pert_inner)
KOOE_outer = column_row_map_space(pert_outer)
KOOE_inner = column_row_map_space(pert_inner)
OE_outer = KOOE_outer & !KO_outer
OE_inner = KOOE_inner & !KO_inner
# make sure OE is kept separate. should be TRUE
any(OE_outer != 0)
any(OE_inner != 0)
write.table(KOOE_outer, file="KOOE_outer.csv", sep=",", quote=F)
write.table(KOOE_inner, file="KOOE_inner.csv", sep=",", quote=F)
# we add NaN entries to J, which will result in their values being ignored in the SSE calculation during gradient descent
J_outer = as.data.table(KOOE_outer, "ORF"); setkey(J_outer, "ORF")
J_inner = as.data.table(KOOE_inner, "ORF"); setkey(J_inner, "ORF")
J_outer[is.na(pert_outer)] = 1
J_inner[is.na(pert_inner)] = 1
# write.table(J_outer, file="J_outer.csv", sep=",", quote=F, row.names=F)
# write.table(J_inner, file="J_inner.csv", sep=",", quote=F, row.names=F)
write.table(J_outer[,!"ORF"], file="J_outer.mat", sep=" ", quote=F, row.names=F, col.names=F)
write.table(J_inner[,!"ORF"], file="J_inner.mat", sep=" ", quote=F, row.names=F, col.names=F)


# make them logical
KO_outer   = KO_outer   == 1
KOOE_outer = KOOE_outer == 1
OE_outer   = OE_outer   == 1
KO_inner   = KO_inner   == 1
KOOE_inner = KOOE_inner == 1
OE_inner   = OE_inner   == 1

# we have to get rid of NaNs. We can try to copy logFC for a node under the same exp conditions
sum(is.na(pert_outer))
pert_outer = as.matrix(pert_outer, "ORF"); pert_outer_updated = pert_outer
pert_inner = as.matrix(pert_inner, "ORF"); pert_inner_updated = pert_inner
pert_outer_updated = NA2avg(pert_outer_updated, pert_inner_updated)
# the remaining NA values are set to zero.
sum(is.na(pert_outer_updated))
sum(is.na(pert_inner_updated))
# correct remaining nodes to zero
pert_outer_updated[is.na(pert_outer_updated)] = 0
pert_inner_updated[is.na(pert_inner_updated)] = 0
# reducing KOs
pert_outer_updated[KO_outer] = pert_outer_updated[KO_outer] - 4
pert_inner_updated[KO_inner] = pert_inner_updated[KO_inner] - 4
# increase OE
pert_outer_updated[OE_outer] = pert_outer_updated[OE_outer] + 1
pert_inner_updated[OE_inner] = pert_inner_updated[OE_inner] + 1
# write. there's spaces in colnames (OE)
# write.table(pert_outer_updated, "logFC_outer.csv", sep=",", quote=F)
# write.table(pert_inner_updated, "logFC_inner.csv", sep=",", quote=F)
write.table(pert_outer_updated, "logFC_outer.mat", sep=" ", quote=F, row.names=F, col.names=F)
write.table(pert_inner_updated, "logFC_inner.mat", sep=" ", quote=F, row.names=F, col.names=F)



## plotting
# plot perturbation values
labels = c("affected", "KO", "OE", "corrected KO", "corrected OE")
color_palette = c("black", "pink", "cyan", "red", "blue")
df1 = data.frame(logFC=pert_outer[!as.matrix(J_outer, "ORF")], label="affected")
df2 = data.frame(logFC=na.omit(pert_outer[KO_outer]), label="KO")
df3 = data.frame(logFC=na.omit(pert_outer[OE_outer]), label="OE")
df4 = data.frame(logFC=pert_outer_updated[KO_outer], label="corrected KO")
df5 = data.frame(logFC=pert_outer_updated[OE_outer], label="corrected OE")
nrow(df1); nrow(df2); nrow(df3); nrow(df4); nrow(df5)
plotdf = rbind(df1, df2, df3, df4, df5)
plotdf$label = factor(plotdf$label, levels=labels)

plt = ggplot(plotdf, aes(logFC, fill=label)) + 
    geom_histogram(alpha=0.5, position="identity", binwidth=0.3)

stepdf = ggplot_build(plt)$data[[1]][,c("xmin", "y", "group")]
stepdf$label = factor(stepdf$group, labels=labels)
# make visual corrections
for (i in 1:length(labels)) {
    limits = c(min(stepdf$xmin[stepdf$group==i])-.1, max(stepdf$xmin[stepdf$group==i])+.1)
    stepdf = rbind(stepdf, data.frame(xmin=limits, y=c(0,0), group=c(i,i), label=labels[i]))
}

plt + geom_step(data=stepdf, aes(x=xmin, y=y, color=label)) +
    theme_linedraw() +
    scale_color_manual(values=color_palette) +
    scale_fill_manual(values=color_palette) +
    scale_x_continuous(name="log fold-change", breaks=c(-10,-5,-1,0,1,5), labels=c("-10","-5","-1","0","1","5")) +
    theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
    scale_y_log10(limits=c(1,1.2e7), expand=c(0,0)) +
    ggtitle("Perturbation corrections") +
    ylab("measurements")

ggsave("perturbation_corrections.pdf", width=7, height=2, units="in")


mydensityplot = function() {
    ggplot(plotdf, aes(logFC, y=..scaled.., fill=label)) + 
        geom_density(alpha=0.7) +
        scale_x_continuous(name="log fold-change", breaks=c(-10,-5,-1,0,1), labels=c("-10","-5","-1","0","1"), limits=c(-12,4)) +
        ylab("scaled density") +
        scale_fill_manual(values=c("pink", "cyan", "red", "blue", "black")) +
        theme_linedraw() +
        theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
        ggtitle("Perturbation corrections")
}

myborderlesshistplot = function() {
    ggplot(plotdf, aes(logFC, fill=label)) + 
        geom_histogram(alpha=0.6, position="identity", bins=80) +
        scale_y_log10(expand = c(0,0)) +
        scale_x_continuous(name="log fold-change", breaks=c(-10,-5,-1,0,1,5), labels=c("-10","-5","-1","0","1")) +
        scale_fill_manual(values=c("black", "pink", "cyan", "red", "blue")) +
        theme_linedraw() +
        theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
        ggtitle("Perturbation corrections")
}

