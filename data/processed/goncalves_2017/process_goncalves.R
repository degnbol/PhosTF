
setwd("~/cwd/data/processed/goncalves_2017")


KOKP = read.table("KP_KO_KP.tsv", row.names=1, sep="\t", header=T, quote="")
KOTF = read.table("KP_KO_TF.tsv", row.names=1, sep="\t", header=T, quote="")
