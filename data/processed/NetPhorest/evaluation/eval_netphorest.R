#!/usr/bin/env Rscript

library(data.table)

flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/processed/NetPhorest/evaluation")

scores = fread("../scores.tsv")
# KP = flatten(fread("~/cwd/data/network/KP.txt", header=F))
# TF = flatten(fread("~/cwd/data/network/TF.txt", header=F))
V = flatten(fread("~/cwd/data/network/V.txt", header=F))
P_eval = fread("~/cwd/data/evaluation/P_eval.tsv")
colnames(P_eval)[colnames(P_eval)=="Source"] = "KP"
P_eval_noknownsite = flatten(fread("~/cwd/data/evaluation/KP_targets_noknownsite.txt", header=F))
KP_edges = fread("~/cwd/data/inference/74/KP_edges.tsv")
protein_KPTF_table = read.table("~/cwd/data/network/node_attributes_full.txt", sep="\t", header=T)
protein_KP = as.character(protein_KPTF_table$ORF[protein_KPTF_table$protein_type0=="KP" & !is.na(protein_KPTF_table$protein_type0)])
protein_TF = as.character(protein_KPTF_table$ORF[protein_KPTF_table$protein_type0=="TF" & !is.na(protein_KPTF_table$protein_type0)])
protein_KPTF = c(protein_KP, protein_TF)
KP = protein_KP
TF = protein_TF

# remove diagonals
P_eval = P_eval[as.character(P_eval$KP) != as.character(P_eval$Target),]
KP_edges = KP_edges[as.character(KP_edges$KP) != as.character(KP_edges$Target),]
# remove non-protein
P_eval = P_eval[as.character(P_eval$Target) %in% protein_KPTF,]
P_eval = P_eval[as.character(P_eval$KP) %in% protein_KPTF,]
KP_edges = KP_edges[as.character(KP_edges$Target) %in% protein_KPTF,]
KP_edges = KP_edges[as.character(KP_edges$KP) %in% protein_KPTF,]



get_evaluation_data = function() {
    venndata = data.frame(
        potential = TRUE,
        biogrid = P_eval$biogrid != "",
        fasolo = !is.na(P_eval$fasolo),
        parca = !is.na(P_eval$parca),
        fiedler = P_eval$fiedler != "",
        yeastkid = !is.na(P_eval$yeastkid) & P_eval$yeastkid > 4.52,
        ptmod = !is.na(P_eval$ptmod) & P_eval$ptmod > 250,
        ptacek = !is.na(P_eval$ptacek)
    )
    venndata$literature = venndata$biogrid | venndata$fasolo | venndata$parca | venndata$fiedler
    venndata$curated = venndata$yeastkid | venndata$ptmod
    venndata$known = venndata$literature | venndata$curated
    venndata$invitro = venndata$known | venndata$ptacek
    venndata$with_site = venndata$invitro & !P_eval$Target%in%P_eval_noknownsite
    venndata
}

venndata = get_evaluation_data()


p.selection = function(col_select, row_idx) {
    drawn = top_marked[row_idx]
    if(missing(row_idx)) row_idx = rep(TRUE, nrow(venndata))
    evaldata = venndata[row_idx, col_select]
    q = sum(evaldata[drawn])  # number of correct inferences
    m = sum(evaldata) # number of true edges
    n = length(evaldata) - m  # number of potential edges (not in true edge set)
    list(value=q, p=phyper(q, m, n, sum(drawn), lower.tail=FALSE))
}



KP2KP.idx = P_eval$Target%in%KP
KP2TF.idx = P_eval$Target%in%TF
stopifnot(all(KP2KP.idx|KP2TF.idx))

# test netphorest drawing similar counts of the KP2KP and KP2TF potential edges  

netphorest_KPs = KP[KP%in%scores$KP]
netphorest_targets = c(KP,TF)[c(KP,TF)%in%scores$Target]
KP2KP_netpho.idx = KP2KP.idx & KP_edges$KP%in%netphorest_KPs & KP_edges$Target%in%netphorest_targets
KP2TF_netpho.idx = KP2TF.idx & KP_edges$KP%in%netphorest_KPs & KP_edges$Target%in%netphorest_targets

ndrawn_KP2KP = sum(KP_edges[KP2KP_netpho.idx, q] < .05)
ndrawn_KP2TF = sum(KP_edges[KP2TF_netpho.idx, q] < .05)
nscored_KP2KP = sum(!is.na(P_eval[KP2KP_netpho.idx,netphorest]))
nscored_KP2TF = sum(!is.na(P_eval[KP2TF_netpho.idx,netphorest]))
drawn_KP2KP = P_eval[KP2KP_netpho.idx,order(-netphorest)[1:ndrawn_KP2KP]]
drawn_KP2TF = P_eval[KP2TF_netpho.idx,order(-netphorest)[1:ndrawn_KP2TF]]
venndata_KP2KP = venndata[KP2KP_netpho.idx, "with_site"]
venndata_KP2TF = venndata[KP2TF_netpho.idx, "with_site"]
q_KP2KP = sum(venndata_KP2KP[drawn_KP2KP])
q_KP2TF = sum(venndata_KP2TF[drawn_KP2TF])
m_KP2KP = sum(venndata_KP2KP)
m_KP2TF = sum(venndata_KP2TF)
n_KP2KP = sum(!venndata_KP2KP)
n_KP2TF = sum(!venndata_KP2TF)

phyper(q_KP2KP, m_KP2KP, n_KP2KP, ndrawn_KP2KP, lower.tail=FALSE)
phyper(q_KP2TF, m_KP2TF, n_KP2TF, ndrawn_KP2TF, lower.tail=FALSE)

# total space of potential edges
n_KP2KP + n_KP2TF + m_KP2KP + m_KP2TF
sum(KP2KP_netpho.idx | KP2TF_netpho.idx)
length(netphorest_KPs) * length(netphorest_targets) # includes self-loops
stopifnot(!any(KP2KP_netpho.idx & KP2TF_netpho.idx))
