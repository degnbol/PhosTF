
# packages
library(reshape2)
library(Matrix)

# functions
flatten = function(x) as.vector(as.matrix(x))
read.vector = function(x) flatten(read.table(x))
read.matrix = function(x) as.matrix(read.table(x))
melt_matrix = function(x) {
    out = melt(as.matrix(x))
    colnames(out) = c("rownames", "colnames", "value")
    out
}

setwd("/Users/christian/GoogleDrev/PKTFX/data/evaluation")

# read
STRING = read.table("../processed/STRING/scores.tsv", sep="\t", header=T, quote="")
STRING_undirected = read.table("../processed/STRING/interactions.tsv", sep="\t", header=T, quote="")
P_biogrid = read.table("../processed/biogrid/P_edges.tsv", sep="\t", header=T, quote="")
P_fasolo = read.table("../processed/fasolo_2011/P_edges_ORF.tsv", sep="\t", header=T, quote="")
P_fiedler = read.table("../processed/fiedler_2009/P_edges_EMAP.tsv", sep="\t", header=T, check.names=F, quote="")[,-c(2,4)]
P_parca = read.table("../processed/parca_2019/P_edges.tsv", sep="\t", header=T, quote="")
P_workman = read.table("../processed/workman_2006/PK_edges.tsv", sep="\t", header=T, quote="")
P_yeastkid = read.table("../processed/yeastkid/P_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,-c(2,4)]
P_netphorest = read.table("../processed/NetPhorest/scores.tsv", sep="\t", header=T, quote="")
P_networkin = read.table("../processed/NetworKIN/scores.tsv", sep="\t", header=T, quote="")[,c(1,2,3,5)] # there is a PCC=.9999 for netphorest so we don't need it
P_networkin_biogrid = read.table("../processed/NetworKIN/scores_biogrid.tsv", sep="\t", header=T, quote="")[,1:3]
T_balaji = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", header=T, quote="")
T_workman = read.table("../processed/workman_2006/TF_edges_nodup.tsv", sep="\t", header=T, quote="")[,c(1,2,4)]
T_yeastract_binding = read.table("../processed/yeastract/binding_ORF.tsv", sep="\t", col.names=c("TF", "Target"), quote="")
T_yeastract_expression = read.table("../processed/yeastract/expression_ORF.tsv", sep="\t", col.names=c("TF", "Target", "Mode"), quote="")
T_harbison = read.table("../processed/harbison_2004/TF_edges.tsv", header=T, sep="\t", quote="")
T_horak = read.table("../processed/horak_2002/TF_edges.tsv", header=T, sep="\t", quote="", check.names=F)[,1:3]
# rename columns to indicate data source or add columns indicating such
colnames(STRING_undirected) = c("Source", "Target", "undirected")
colnames(P_biogrid) = c("Source", "Target", "biogrid")
P_fasolo$fasolo = 1; colnames(P_fasolo) = c("Source", "Target", "fasolo")
colnames(P_fiedler) = c("Source", "Target", "fiedler", "EMAP")
P_parca$parca = 1; colnames(P_parca) = c("Source", "Target", "parca")
colnames(P_workman) = c("Source", "Target", "workman SLL", "workman pval")
colnames(P_yeastkid) = c("Source", "Target", "yeastkid")
colnames(P_netphorest) = c("Source", "Target", "netphorest")
colnames(P_networkin) = c("Source", "Target", "networkin", "networkin STRING")
colnames(P_networkin_biogrid) = c("Source", "Target", "networkin_biogrid")
T_balaji$balaji = 1; colnames(T_balaji) = c("Source", "Target", "balaji")
colnames(T_workman) = c("Source", "Target", "workman")
T_yeastract_binding$"yeastract binding" = 1; colnames(T_yeastract_binding) = c("Source", "Target", "yeastract binding")
colnames(T_yeastract_expression) = c("Source", "Target", "yeastract expression")
colnames(T_harbison) = c("Source", "Target", "harbison")
colnames(T_horak) = c("Source", "Target", "horak")

KP_fname = "../perturbation/KP.txt"
TF_fname = "../perturbation/TF.txt"
PTX_fname = "../perturbation/PTX.txt"
KP = read.vector(KP_fname)
TF = read.vector(TF_fname)
PTX = read.vector(PTX_fname)
PT = c(KP,TF)

# make source and target identification that will match the melted version of inferred adjacencny matrices
P_eval = melt_matrix(sparseMatrix(c(), c(), dims=c(length(PT),  length(KP)), dimnames=list(PT, KP)))[,c(2,1)]
T_eval = melt_matrix(sparseMatrix(c(), c(), dims=c(length(PTX), length(TF)), dimnames=list(PTX, TF)))[,c(2,1)]
colnames(P_eval) = c("Source", "Target")
colnames(T_eval) = c("Source", "Target")

# merge tables
P_eval = merge(P_eval, P_biogrid, all.x=T)
P_eval = merge(P_eval, P_fasolo, all.x=T)
P_eval = merge(P_eval, P_parca, all.x=T)
P_eval = merge(P_eval, P_fiedler, all.x=T)
P_eval = merge(P_eval, P_workman, all.x=T)
P_eval = merge(P_eval, P_yeastkid, all.x=T)
P_eval = merge(P_eval, P_netphorest, all.x=T)
P_eval = merge(P_eval, P_networkin, all.x=T)
P_eval = merge(P_eval, P_networkin_biogrid, all.x=T)
T_eval = merge(T_eval, T_balaji, all.x=T)
T_eval = merge(T_eval, T_workman, all.x=T)
T_eval = merge(T_eval, T_yeastract_binding, all.x=T)
T_eval = merge(T_eval, T_yeastract_expression, all.x=T)
T_eval = merge(T_eval, T_yeastract_expression, all.x=T)
T_eval = merge(T_eval, T_harbison, all.x=T)
T_eval = merge(T_eval, T_horak, all.x=T)

P_eval = merge(P_eval, STRING, all.x=T)
T_eval = merge(T_eval, STRING, all.x=T)
P_eval = merge(P_eval, STRING_undirected, all.x=T)
T_eval = merge(T_eval, STRING_undirected, all.x=T)

P_eval$n_datasets = rowSums(!is.na(P_eval[,colnames(P_eval) != "workman pval"][,-c(1,2)]))
T_eval$n_datasets = rowSums(!is.na(T_eval[,colnames(T_eval) != "workman pval"][,-c(1,2)]))
sum(P_eval$n_datasets > 0)
sum(T_eval$n_datasets > 0)
# number of edges with directed info
sum(rowSums(!is.na(P_eval[,!(colnames(P_eval) %in% c("workman pval", "undirected", "n_datasets"))][,-c(1,2)])) > 0)
sum(rowSums(!is.na(T_eval[,!(colnames(T_eval) %in% c("workman pval", "undirected", "n_datasets"))][,-c(1,2)])) > 0)


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

positives1 = (P_eval$yeastkid > 4.52) | (P_eval$undirected > 950) | (P_eval$ptmod > mean(P_eval$ptmod, na.rm=T)) | P_eval$fasolo | P_eval$parca | !is.na(P_eval$fiedler) | !is.na(P_eval$biogrid) | P_eval$`workman pval` < .05 | P_eval$netphorest > .27 | P_eval$networkin_biogrid > 20
negatives1 = (P_eval$yeastkid < 0.00) | (P_eval$undirected < 155) | (P_eval$ptmod < mean(P_eval$ptmod, na.rm=T)) | P_eval$netphorest < .04 | P_eval$networkin_biogrid < .0002 
positives2 = (P_eval$yeastkid > 4.52) | (P_eval$ptmod > mean(P_eval$ptmod, na.rm=T)) | P_eval$fasolo | P_eval$parca | !is.na(P_eval$fiedler) | !is.na(P_eval$biogrid) | P_eval$`workman pval` < .05 
negatives2 = (P_eval$yeastkid < 0.00) | (P_eval$ptmod < mean(P_eval$ptmod, na.rm=T))
positives3 = P_eval$ptmod > mean(P_eval$ptmod, na.rm=T) | P_eval$n_datasets >= 7
negatives3 = P_eval$ptmod < mean(P_eval$ptmod, na.rm=T)
qntl = .9
positives4 = 
    (rowSums(!is.na(P_eval[,c("fasolo", "parca", "fiedler", "biogrid")])) > 1) |
    (P_eval$yeastkid > 4.52 | is.na(P_eval$yeastkid)) &
    (P_eval$undirected > quantile(P_eval$undirected, qntl, na.rm=T) | is.na(P_eval$undirected)) &
    (P_eval$ptmod > mean(P_eval$ptmod, na.rm=T) | is.na(P_eval$ptmod)) &
    (P_eval$networkin_biogrid > quantile(P_eval$networkin_biogrid, qntl, na.rm=T) | is.na(P_eval$networkin_biogrid)) &
    (P_eval$netphorest > quantile(P_eval$netphorest, qntl, na.rm=T) | is.na(P_eval$netphorest)) &
    (P_eval$n_datasets > 0) # otherwise the initial logic holds true for edges where everything is NA
negatives4 = 
    is.na(P_eval$fasolo) &
    is.na(P_eval$parca) &
    is.na(P_eval$fiedler) &
    is.na(P_eval$biogrid) &
    (P_eval$yeastkid < mean(P_eval$yeastkid, na.rm=T) | is.na(P_eval$yeastkid)) &
    (P_eval$undirected < quantile(P_eval$undirected, 1-qntl, na.rm=T) | is.na(P_eval$undirected)) &
    (P_eval$ptmod < mean(P_eval$ptmod, na.rm=T) | is.na(P_eval$ptmod)) &
    (P_eval$networkin_biogrid < quantile(P_eval$networkin_biogrid, 1-qntl, na.rm=T) | is.na(P_eval$networkin_biogrid)) &
    (P_eval$netphorest < quantile(P_eval$netphorest, 1-qntl, na.rm=T) | is.na(P_eval$netphorest)) &
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
write.table(T_eval, file="T_eval.tsv", sep="\t", quote=F, row.names=F, na="")



