
library(stringr)

# functions
flatten = function(x) as.vector(as.matrix(x))

# get values where column name == row name
column_row_values = function(mat) {
    long = melt(mat)
    long$value[as.character(long$Var1) == as.character(long$Var2)]
}

setwd("~/cwd/data/processed/chua_2006")
# check names off so duplicate names are not changed
TFOE = as.matrix(read.table("OE.tsv", header=T, row.names=1, check.names=F))
TFKO = as.matrix(read.table("KO.tsv", header=T, row.names=1, check.names=F))

TFOE_plus = TFOE[,str_sub(colnames(TFOE), start=-1) == "+"]; colnames(TFOE_plus) = gsub("\\+$", "", colnames(TFOE_plus))
TFOE_minu = TFOE[,str_sub(colnames(TFOE), start=-1) == "-"]; colnames(TFOE_minu) = gsub("\\-$", "", colnames(TFOE_minu))
TFKO_plus = TFKO[,str_sub(colnames(TFKO), start=-1) == "+"]; colnames(TFKO_plus) = gsub("\\+$", "", colnames(TFKO_plus))
TFKO_minu = TFKO[,str_sub(colnames(TFKO), start=-1) == "-"]; colnames(TFKO_minu) = gsub("\\-$", "", colnames(TFKO_minu))
TFOE_na = is.na(TFOE_plus) | is.na(TFOE_minu)
TFKO_na = is.na(TFKO_plus) | is.na(TFKO_minu)
TFOE_plus_nona = TFOE_plus[!TFOE_na]
TFOE_minu_nona = TFOE_minu[!TFOE_na]
TFKO_plus_nona = TFKO_plus[!TFKO_na]
TFKO_minu_nona = TFKO_minu[!TFKO_na]
cor(TFOE_plus_nona, TFOE_minu_nona)
cor(TFKO_plus_nona, TFKO_minu_nona)
# seems we do not need to reverse sign and seems that ko has low correlation between the two replicates

# look at values for the pertubated nodes
OE_plus = column_row_values(TFOE_plus)
OE_minu = column_row_values(TFOE_minu)
KO_plus = column_row_values(TFKO_plus)
KO_minu = column_row_values(TFKO_minu)

plot(density(OE_plus[!is.na(OE_plus)]))
plot(density(OE_minu[!is.na(OE_minu)]))
plot(density(KO_plus[!is.na(KO_plus)]))
plot(density(KO_minu[!is.na(KO_minu)]))
hist(OE_plus, breaks=10)
hist(OE_minu, breaks=10)
hist(KO_plus, breaks=10)
hist(KO_minu, breaks=10)

OE_na = is.na(OE_plus) | is.na(OE_minu)
KO_na = is.na(KO_plus) | is.na(KO_minu)
cor(OE_plus[!OE_na], OE_minu[!OE_na])
cor(KO_plus[!KO_na], KO_minu[!KO_na])


# KO has high correlation between replicates for the KO gene, which is good and expected.
# so the data is clearly signed correctly already











