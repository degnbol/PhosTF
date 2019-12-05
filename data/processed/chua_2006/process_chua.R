# functions
flatten = function(x) as.vector(as.matrix(x))

# get values where column name == row name
column_row_values = function(x) {
    values = c()
    for (j in 1:ncol(x)) {
        values = c(values, x[colnames(x)[j] == rownames(x),j])
    }
    values
}

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
TFOE1_nona = TFOE1[!TFOE_na]
TFOE2_nona = TFOE2[!TFOE_na]
TFKO1_nona = TFKO1[!TFKO_na]
TFKO2_nona = TFKO2[!TFKO_na]

cor(TFOE1_nona, TFOE2_nona)
cor(TFKO1_nona, TFKO2_nona)
# seems we do not need to reverse sign and seems that ko has low correlation between the two replicates

# look at values for the pertubated nodes
OE1 = column_row_values(TFOE1)
OE2 = column_row_values(TFOE2)
KO1 = column_row_values(TFKO1)
KO2 = column_row_values(TFKO2)

plot(density(OE1[!is.na(OE1)]))
plot(density(OE2[!is.na(OE2)]))
plot(density(KO1[!is.na(KO1)]))
plot(density(KO2[!is.na(KO2)]))
hist(OE1, breaks=10)
hist(OE2, breaks=10)
hist(KO1, breaks=10)
hist(KO2, breaks=10)

OE_na = is.na(OE1) | is.na(OE2)
KO_na = is.na(KO1) | is.na(KO2)
cor(OE1[!OE_na], OE2[!OE_na])
cor(KO1[!KO_na], KO2[!KO_na])


# KO has high correlation between replicates for the KO gene, which is good and expected.
# so the data is clearly signed correctly already











