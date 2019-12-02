# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("/Users/christian/GoogleDrev/PKTFX/data/processed/chua_2006")
# check names off so duplicate names are not changed
TFOE = read.table("TF_OE.tsv", header=T, row.names=1, check.names=F)
TFKO = read.table("TF_KO.tsv", header=T, row.names=1, check.names=F)

TFOE1 = TFOE[,seq(1,ncol(TFOE),2)]
TFOE2 = TFOE[,seq(2,ncol(TFOE),2)]
TFKO1 = TFKO[,seq(1,ncol(TFKO),2)]
TFKO2 = TFKO[,seq(2,ncol(TFKO),2)]
TFOE_na = is.na(TFOE1) | is.na(TFOE2)
TFKO_na = is.na(TFKO1) | is.na(TFKO2)
TFOE1 = TFOE1[!TFOE_na]
TFOE2 = TFOE2[!TFOE_na]
TFKO1 = TFKO1[!TFKO_na]
TFKO2 = TFKO2[!TFKO_na]

cor(TFOE1, TFOE2)
cor(TFKO1, TFKO2)
# seems we do not need to reverse sign and seems that ko has low correlation between the two replicates

# look at values for the pertubated nodes
OE = rep(0, ncol(TFOE))
KO = rep(0, ncol(TFKO))
OEidx = match(colnames(TFOE), rownames(TFOE))
KOidx = match(colnames(TFKO), rownames(TFKO))
for (i in 1:ncol(TFKO)) {KO[i] = TFKO[rowidx[i],i]}
for (i in 1:ncol(TFOE)) {OE[i] = TFOE[rowidx[i],i]}
KO1 = KO[seq(1,length(KO),2)]
KO2 = KO[seq(2,length(KO),2)]
OE1 = OE[seq(1,length(OE),2)]
OE2 = OE[seq(2,length(OE),2)]
OE_na = is.na(OE1) | is.na(OE2)
KO_na = is.na(KO1) | is.na(KO2)
cor(OE1[!OE_na], OE2[!OE_na])
cor(KO1[!KO_na], KO2[!KO_na])
hist(OE, breaks=10)
hist(KO, breaks=10)
hist(OE1, breaks=10)
hist(OE2, breaks=10)
hist(KO1, breaks=10)
hist(KO2, breaks=10)

# KO has high correlation between replicates for the KO gene, which is good and expected.
# so the data is clearly signed correctly already






