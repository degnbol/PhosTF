#!/usr/bin/env Rscript
library(data.table)
library(reshape2)
library(Matrix)
suppressPackageStartupMessages(library(here))
ROOT = paste0(here(), "/")
cwd = function(path) paste0(ROOT, path)

# functions
flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
read.matrix = function(x) as.matrix(read.table(x))
melt_matrix = function(x) {
    out = melt(as.matrix(x))
    colnames(out) = c("rownames", "colnames", "value")
    out
}

# read
STRING = read.table(cwd("data/STRING/scores.tsv"), sep="\t", header=T, quote="")
STRING_undirected = read.table(cwd("data/STRING/interactions.tsv"), sep="\t", header=T, quote="")
P_biogrid = read.table(cwd("data/biogrid/P_edges.tsv"), sep="\t", header=T, quote="")
P_fasolo = read.table(cwd("data/fasolo_2011/P_edges_ORF.tsv"), sep="\t", header=T, quote="")
P_fiedler = read.table(cwd("data/fiedler_2009/P_edges_EMAP.tsv"), sep="\t", header=T, check.names=F, quote="")[,-c(2,4)]
P_parca = read.table(cwd("data/parca_2019/P_edges.tsv"), sep="\t", header=T, quote="")
# P_workman = read.table(cwd("data/workman_2006/PK_edges.tsv"), sep="\t", header=T, quote="")
P_yeastkid = read.table(cwd("data/yeastkid/P_edges.tsv"), sep="\t", header=T, quote="", check.names=F)[,-c(2,4)]
P_netphorest = read.table(cwd("data/NetPhorest/scores.tsv"), sep="\t", header=T, quote="")
P_networkin = read.table(cwd("data/NetworKIN/scores.tsv"), sep="\t", header=T, quote="")[,c(1,2,3,5)] # there is a PCC=.9999 for netphorest so we don't need it
P_networkin_biogrid = read.table(cwd("data/NetworKIN/scores_biogrid.tsv"), sep="\t", header=T, quote="")[,1:3]
P_ptacek = read.table(cwd("data/ptacek_2005/KP_edges.tsv"), sep="\t", header=T, quote="", check.names=F)
# T_balaji = read.table(cwd("data/balaji_2006/TF_edges.tsv"), sep="\t", header=T, quote="")
# T_workman = read.table(cwd("data/workman_2006/TF_edges.tsv"), sep="\t", header=T, quote="")[,c(1,2,4)]
# T_yeastract_binding = read.table(cwd("data/yeastract/binding_ORF.tsv"), sep="\t", col.names=c("TF", "Target"), quote="")
# T_yeastract_expression = read.table(cwd("data/yeastract/expression_ORF.tsv"), sep="\t", col.names=c("TF", "Target", "Mode"), quote="")
# T_harbison_YPD = read.table(cwd("data/harbison_2004/TF_edges_YPD.tsv"), header=T, sep="\t", quote="")
# T_harbison_conds = read.table(cwd("data/harbison_2004/TF_edges_conds.tsv"), header=T, sep="\t", quote="")
# T_horak = read.table(cwd("data/horak_2002/TF_edges.tsv"), header=T, sep="\t", quote="", check.names=F)[,1:3]
# rename columns to indicate data source or add columns indicating such
colnames(STRING_undirected) = c("Source", "Target", "undirected")
colnames(P_biogrid) = c("Source", "Target", "biogrid")
P_fasolo$fasolo = 1; colnames(P_fasolo) = c("Source", "Target", "fasolo")
colnames(P_fiedler) = c("Source", "Target", "fiedler", "EMAP")
colnames(P_parca) = c("Source", "Target"); P_parca$parca = 1
# colnames(P_workman) = c("Source", "Target", "workman SLL", "workman pval")
colnames(P_yeastkid) = c("Source", "Target", "yeastkid")
colnames(P_netphorest) = c("Source", "Target", "netphorest")
colnames(P_networkin) = c("Source", "Target", "networkin", "networkin STRING")
colnames(P_networkin_biogrid) = c("Source", "Target", "networkin_biogrid")
colnames(P_ptacek) = c("Source", "Target"); P_ptacek$ptacek = 1
# colnames(T_balaji) = c("Source", "Target", "balaji")
# colnames(T_workman) = c("Source", "Target", "workman")
# T_yeastract_binding$"yeastract binding" = 1; colnames(T_yeastract_binding) = c("Source", "Target", "yeastract binding")
# colnames(T_yeastract_expression) = c("Source", "Target", "yeastract expression")
# colnames(T_harbison_YPD) = c("Source", "Target", "harbison_YPD")
# colnames(T_harbison_conds) = c("Source", "Target", "harbison_conds")
# colnames(T_horak) = c("Source", "Target", "horak")

KP = fread(cwd("KP_protein.tsv"))$ORF
TF = read.vector(cwd("TF.txt"))
V = read.vector(cwd("V_protein.txt"))
PT = c(KP,TF)

# make source and target identification that will match the melted version of inferred adjacencny matrices
P_eval = melt_matrix(sparseMatrix(c(), c(), dims=c(length(PT), length(KP)), dimnames=list(PT, KP)))[,c(2,1)]
# T_eval = melt_matrix(sparseMatrix(c(), c(), dims=c(length(V),  length(TF)), dimnames=list(V,  TF)))[,c(2,1)]
colnames(P_eval) = c("Source", "Target")
# colnames(T_eval) = c("Source", "Target")

# merge tables
P_eval = merge(P_eval, P_biogrid, all.x=T)
P_eval = merge(P_eval, P_fasolo, all.x=T)
P_eval = merge(P_eval, P_parca, all.x=T)
P_eval = merge(P_eval, P_fiedler, all.x=T)
# P_eval = merge(P_eval, P_workman, all.x=T)
P_eval = merge(P_eval, P_yeastkid, all.x=T)
P_eval = merge(P_eval, P_netphorest, all.x=T)
P_eval = merge(P_eval, P_networkin, all.x=T)
P_eval = merge(P_eval, P_networkin_biogrid, all.x=T)
P_eval = merge(P_eval, STRING, all.x=T)
P_eval = merge(P_eval, STRING_undirected, all.x=T)
P_eval = merge(P_eval, P_ptacek, all.x=T)

# T_eval = merge(T_eval, T_balaji, all.x=T)
# T_eval = merge(T_eval, T_workman, all.x=T)
# T_eval = merge(T_eval, T_yeastract_binding, all.x=T)
# T_eval = merge(T_eval, T_yeastract_expression, all.x=T)
# T_eval = merge(T_eval, T_harbison_YPD, all.x=T)
# T_eval = merge(T_eval, T_harbison_conds, all.x=T)
# T_eval = merge(T_eval, T_horak, all.x=T)
# T_eval = merge(T_eval, STRING, all.x=T)
# T_eval = merge(T_eval, STRING_undirected, all.x=T)

P_eval$n_datasets = rowSums(!is.na(P_eval[,-c(1,2)]))
# T_eval$n_datasets = rowSums(!is.na(T_eval[,-c(1,2)]))



sum(!is.na(P_eval$biogrid))
sum(!is.na(P_eval$fasolo))
sum(!is.na(P_eval$parca))
sum(!is.na(P_eval$fiedler))
sum(!is.na(P_eval$ptacek))
sum(P_eval$n_datasets > 0)
# sum(T_eval$n_datasets > 0)
# number of edges with directed info
sum(rowSums(!is.na(P_eval[,!(colnames(P_eval) %in% c("undirected", "n_datasets"))][,-c(1,2)])) > 0)
# sum(rowSums(!is.na(T_eval[,!(colnames(T_eval) %in% c("undirected", "n_datasets"))][,-c(1,2)])) > 0)

# stratification
# summary(P_eval$undirected, na.rm=T)
# summary(P_eval$yeastkid, na.rm=T)
# summary(P_eval$ptmod, na.rm=T)
# summary(P_eval$expression, na.rm=T)

# sum(P_eval$yeastkid > 4.52, na.rm=T)
# sum(P_eval$yeastkid < 0, na.rm=T)
# sum(P_eval$undirected > 950, na.rm=T)
# sum(P_eval$undirected < 155, na.rm=T)
# sum(P_eval$ptmod > 300, na.rm=T)
# sum(P_eval$ptmod < 300, na.rm=T)

ptmod_pos = P_eval$ptmod > mean(P_eval$ptmod, na.rm=T)
ptmod_neg = P_eval$ptmod < mean(P_eval$ptmod, na.rm=T)
yeastkid_pos = P_eval$yeastkid > 4.52
yeastkid_neg = P_eval$yeastkid < 0.00
positives1 = yeastkid_pos | (P_eval$undirected > 950) | ptmod_pos | P_eval$fasolo | P_eval$parca | !is.na(P_eval$fiedler) | !is.na(P_eval$biogrid) | !is.na(P_eval$ptacek)
negatives1 = yeastkid_neg | (P_eval$undirected < 155) | ptmod_neg
positives2 = yeastkid_pos | ptmod_pos | P_eval$fasolo | P_eval$parca | !is.na(P_eval$fiedler) | !is.na(P_eval$biogrid)
negatives2 = yeastkid_neg | ptmod_neg
positives3 = ptmod_pos | P_eval$n_datasets >= 7
negatives3 = ptmod_neg
qntl = .9
positives4 = 
    (rowSums(!is.na(P_eval[,c("fasolo", "parca", "fiedler", "biogrid")])) > 1) |
    (yeastkid_pos | is.na(P_eval$yeastkid)) &
    (P_eval$undirected > quantile(P_eval$undirected, qntl, na.rm=T) | is.na(P_eval$undirected)) &
    (ptmod_pos | is.na(P_eval$ptmod)) &
    (P_eval$networkin_biogrid > quantile(P_eval$networkin_biogrid, qntl, na.rm=T) | is.na(P_eval$networkin_biogrid)) &
    (P_eval$n_datasets > 0) # otherwise the initial logic holds true for edges where everything is NA
negatives4 = 
    is.na(P_eval$fasolo) &
    is.na(P_eval$parca) &
    is.na(P_eval$fiedler) &
    is.na(P_eval$biogrid) &
    (yeastkid_neg | is.na(P_eval$yeastkid)) &
    (P_eval$undirected < quantile(P_eval$undirected, 1-qntl, na.rm=T) | is.na(P_eval$undirected)) &
    (ptmod_neg | is.na(P_eval$ptmod)) &
    (P_eval$networkin_biogrid < quantile(P_eval$networkin_biogrid, 1-qntl, na.rm=T) | is.na(P_eval$networkin_biogrid)) &
    (P_eval$n_datasets > 0) # otherwise the initial logic holds true for edges where everything is NA

# NA means the entry is not in the set
positives1[is.na(positives1)] = F
negatives1[is.na(negatives1)] = F
positives2[is.na(positives2)] = F
negatives2[is.na(negatives2)] = F
positives3[is.na(positives3)] = F
negatives3[is.na(negatives3)] = F
positives4[is.na(positives4)] = F
negatives4[is.na(negatives4)] = F
# set them to zeros and ones
P_eval$goldstandard1 = NA
P_eval$goldstandard2 = NA
P_eval$goldstandard3 = NA
P_eval$goldstandard4 = NA
P_eval$goldstandard1[positives1 & !negatives1] = 1
P_eval$goldstandard1[negatives1 & !positives1] = 0
P_eval$goldstandard2[positives2 & !negatives2] = 1
P_eval$goldstandard2[negatives2 & !positives2] = 0
P_eval$goldstandard3[positives3 & !negatives3] = 1
P_eval$goldstandard3[negatives3] = 0
P_eval$goldstandard4[positives4 & !negatives4] = 1
P_eval$goldstandard4[negatives4 & !positives4] = 0


write.table(P_eval, file="P_eval.tsv", sep="\t", quote=F, row.names=F, na="")
# write.table(T_eval, file="T_eval.tsv", sep="\t", quote=F, row.names=F, na="")



