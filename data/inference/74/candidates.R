#!/usr/bin/env Rscript

"
take KPs from KP->TF with high weight and many shared GO terms
that are actually protein kinase/phosphatase 
hopefully some will be known edges, 
talk about the thing the TF is known to do
many KP->TF edges are negative, hopefully this can be understood to make biological sense in these few examples. 
"



library(data.table)
library(ggplot2)


setwd("~/cwd/data/inference/74")


KP_edges = fread("KP_edges_eval_noTN.tsv", header=T)
shared_GO = fread("~/cwd/data/go/shared_GO.tsv", header=T)

shared_P = shared_GO[Aspect=="P" & Substrate=="TF",]
setnames(shared_P, "ORF", "Target")

edges = na.omit(KP_edges[shared_P, on=c("KP","Target")])[order(-shared, -abs(marker)),]
edges[,Candidate:=shared*abs(marker)]
edges = edges[shared>0,]
edges[,target_type:=NULL]
edges[,Aspect:=NULL]
edges[,Substrate:=NULL]
setnames(edges, "Target", "TF")
edges = edges[order(-Candidate),]

thres = 6
(edges[evaluation_set==T & shared>=thres,.N] / edges[shared>=thres,.N]) / (edges[evaluation_set==T & shared<thres,.N] / edges[shared<thres,.N])

# KP repression or activation (repression since largest abs(marker) are negative)
View(edges[shared>=6,][order(-abs(marker)),])
# number 1 is YER133W	YDR207C which seems to both have something to do with meiotic regulation
View(edges[shared>=5,][order(-abs(marker)),])
# number 2 is YHL007C	YHR084W which are called STE20 and STE12 and are obviously related (pseudohyphal/invasive growth)
View(edges[shared>=4,][order(-abs(marker)),])
# number 1 is YLR362W	YHR084W which are called STE11 and STE12, same stuff as before
View(edges[shared>=3,][order(-abs(marker)),])
# number 1 is still YLR362W	YHR084W
View(edges[shared>=2,][order(-abs(marker)),])
# number 1 is YNR031C	YHR206W indicated as suppressors, which corresponds to the sign of edge weight, both related to osmotic stress
View(edges[shared>=1,][order(-abs(marker)),])
# number 1 again is YNR031C	YHR206W

# KP activation
View(edges[shared>=6,][order(-marker),])
# number 1 is YBR160W	YLR182W. G1 and cell cycle both mentioned to elicit activating roles
View(edges[shared>=5,][order(-marker),])
# number 1 is YBR160W	YKL112W. YKL112W involved with DNA replication
View(edges[shared>=4,][order(-marker),])
# number 1 is YAR019C	YJR094C. YAR019C phosphorylates pol II in meiosis, YJR094C is master regulator of meiosis, activates.
View(edges[shared>=3,][order(-marker),])
# number 1 is YAR019C	YJR094C again
View(edges[shared>=2,][order(-marker),])
# number 1 is YAR019C	YJR094C again
View(edges[shared>=1,][order(-marker),])
# number 1 is YAR019C	YJR094C again