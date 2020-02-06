
# purpose is to make column EDGES in table for edge data, 
# indicating how many unique edges are described in the data based on having source node among KP or TF and target node among V

options(stringsAsFactors=F)

# functions

flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/network")

KPs = flatten(read.table("KP.txt"))
TFs = flatten(read.table("TF.txt"))
Vs = flatten(read.table("V.txt"))
KPTFs = c(KPs,TFs)

biogrid = read.table("../processed/biogrid/P_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
fasolo = read.table("../processed/fasolo_2011/P_edges_ORF.tsv", sep="\t", header=T, quote="", check.names=F)
parca = read.table("../processed/parca_2019/P_edges.tsv", sep="\t", header=T, quote="", check.names=F)
fiedler = read.table("../processed/fiedler_2009/P_edges_ORF.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
fiedler_EMAP = read.table("../processed/fiedler_2009/EMAP_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,c(1,3)]
yeastkid = read.table("../processed/yeastkid/P_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,c(1,3)]
colnames(fiedler) = colnames(fiedler_EMAP) = colnames(yeastkid) = c("P", "Target")
netphorest = read.table("../processed/NetPhorest/scores.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
networkin = read.table("../processed/NetworKIN/scores.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
ptacek = read.table("../processed/ptacek_2005/KP_edges.tsv", sep="\t", header=T, quote="", check.names=F)
balaji = read.table("../processed/balaji_2006/TF_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
workman = read.table("../processed/workman_2006/TF_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
harbison_YPD = read.table("../processed/harbison_2004/TF_edges_YPD.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
harbison_conds = read.table("../processed/harbison_2004/TF_edges_conds.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
horak = read.table("../processed/horak_2002/TF_edges.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
yeastract_binding = read.table("../processed/yeastract/binding_ORF.tsv", sep="\t", quote="", check.names=F, col.names=c("TF", "Target"))
yeastract_activators = read.table("../processed/yeastract/activator_expression_ORF.tsv", sep="\t", quote="", check.names=F, col.names=c("TF", "Target"))
yeastract_inhibitors = read.table("../processed/yeastract/inhibitor_expression_ORF.tsv", sep="\t", quote="", check.names=F, col.names=c("TF", "Target"))
STRING_directed = read.table("../processed/STRING/directed.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]
STRING_undirected = read.table("../processed/STRING/interactions.tsv", sep="\t", header=T, quote="", check.names=F)[,1:2]

yeastract_ambiguous = yeastract_activators[paste(yeastract_activators$TF,yeastract_activators$Target) %in% paste(yeastract_inhibitors$TF,yeastract_inhibitors$Target),]
yeastract_expression = rbind(yeastract_activators, yeastract_inhibitors)
yeastract_expression = yeastract_expression[!(paste(yeastract_expression$TF, yeastract_expression$Target) %in% paste(yeastract_ambiguous$TF, yeastract_ambiguous$Target)),]

nrow(harbison_YPD)+nrow(harbison_conds)
# we could count unique rows since we only want to count the binding scores from STRING, not catalysis etc. but idk
nrow(unique(STRING_directed))

# make sure not to count self edges or the same edges multiple times
biogrid = unique(biogrid[biogrid$P!=biogrid$Target,])
fasolo = unique(fasolo[fasolo$P!=fasolo$Target,])
parca = unique(parca[parca$P!=parca$Target,])
fiedler = unique(fiedler[fiedler$P!=fiedler$Target,])
fiedler_EMAP = unique(fiedler_EMAP[fiedler_EMAP$P!=fiedler_EMAP$Target,])
yeastkid = unique(yeastkid[yeastkid$P!=yeastkid$Target,])
netphorest = unique(netphorest[netphorest$KP!=netphorest$Target,])
networkin = unique(networkin[networkin$KP!=networkin$Target,])
ptacek = unique(ptacek[ptacek$KP!=ptacek$Target,])
balaji = unique(balaji[balaji$TF!=balaji$Target,])
workman = unique(workman[workman$TF!=workman$Target,])
harbison = unique(rbind(harbison_YPD,harbison_conds))
harbison = harbison[harbison$TF!=harbison$Target,]
horak = unique(horak[horak$TF!=horak$Target,])
yeastract_binding = unique(yeastract_binding[yeastract_binding$TF!=yeastract_binding$Target,])
yeastract_expression = unique(yeastract_expression[yeastract_expression$TF!=yeastract_expression$Target,])
yeastract_ambiguous = unique(yeastract_ambiguous[yeastract_ambiguous$TF!=yeastract_ambiguous$Target,])
STRING_directed = unique(STRING_directed[STRING_directed$Source!=STRING_directed$Target,])
STRING_undirected = unique(STRING_undirected[STRING_undirected$protein1!=STRING_undirected$protein2,])

# counts of edges in each dataset
sum((biogrid$P %in% KPs) & (biogrid$Target %in% KPTFs))
sum((fasolo$P %in% KPs) & (fasolo$Target %in% KPTFs))
sum((parca$P %in% KPs) & (parca$Target %in% KPTFs))
sum((fiedler$P %in% KPs) & (fiedler$Target %in% KPTFs))
sum((fiedler_EMAP$P %in% KPs) & (fiedler_EMAP$Target %in% KPTFs))
sum((yeastkid$P %in% KPs) & (yeastkid$Target %in% KPTFs))
sum((netphorest$KP %in% KPs) & (netphorest$Target %in% KPTFs))
sum((networkin$KP %in% KPs) & (networkin$Target %in% KPTFs))
sum((ptacek$KP %in% KPs) & (ptacek$Target %in% KPTFs))
sum((balaji$TF %in% TFs) & (balaji$Target %in% Vs))
sum((workman$TF %in% TFs) & (workman$Target %in% Vs))
sum((harbison$TF %in% TFs) & (harbison$Target %in% Vs))
sum((horak$TF %in% TFs) & (horak$Target %in% Vs))
sum((yeastract_binding$TF %in% TFs) & (yeastract_binding$Target %in% Vs))
sum((yeastract_expression$TF %in% TFs) & (yeastract_expression$Target %in% Vs))
sum((yeastract_ambiguous$TF %in% TFs) & (yeastract_ambiguous$Target %in% Vs))
sum((STRING_directed$Source %in% KPs) & (STRING_directed$Target %in% KPTFs))
sum((STRING_directed$Source %in% TFs) & (STRING_directed$Target %in% Vs))
sum((STRING_undirected$protein1 %in% TFs) & (STRING_undirected$protein2 %in% Vs))

colnames(netphorest)[1] = "P"
colnames(networkin)[1] = "P"

# counts of total described KP edges versus total potential
KP_edges = unique(rbind(biogrid, fasolo, parca, fiedler, fiedler_EMAP, yeastkid))
KP_edges_pred = unique(rbind(KP_edges, netphorest, networkin, data.frame(P=STRING_directed$Source, Target=STRING_directed$Target)))
colnames(STRING_directed)[1] = "TF"
colnames(STRING_undirected) = c("TF", "Target")
TF_edges = unique(rbind(balaji, harbison, horak, yeastract_binding, yeastract_expression, yeastract_ambiguous, STRING_directed, STRING_undirected))
sum((KP_edges$P %in% KPs) & (KP_edges$Target %in% KPTFs))
sum((KP_edges_pred$P %in% KPs) & (KP_edges_pred$Target %in% KPTFs))
sum((TF_edges$TF %in% TFs) & (TF_edges$Target %in% Vs))
length(KPs) * length(KPTFs) - length(KPs)  # potential (no self-loops)
length(TFs) * length(Vs) - length(TFs)  # potential (no self-loops)






