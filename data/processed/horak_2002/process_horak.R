
setwd("~/cwd/data/processed/horak_2002")

raw = read.table("raw_ORF.tsv", sep="\t", header=T, quote="")
raw$Target = gsub("-([0-9])[A-Z]$", "-\\1", gsub("-[0-9]$", "", raw$Target))

agg = aggregate(Score ~ TF + Target, data=raw, max)
agg$Pval = raw$Pval[match(paste(agg$TF, agg$Target), paste(raw$TF, raw$Target))]
agg$Pval_est = (1-pnorm(agg$Score))*2

# check that the pval estimate makes sense
correct = 
    (agg$Pval == ">0.01" & agg$Pval_est >= 0.01) | 
    (agg$Pval == "<0.01" & agg$Pval_est <= 0.01 & agg$Pval_est >= 0.001) |
    (agg$Pval == "<0.001" & agg$Pval_est <= 0.001 & agg$Pval_est >= 0.0001) |
    (agg$Pval == "<0.0001" & agg$Pval_est <= 0.0001)
sum(correct)
sum(!correct)
# it is the most accurate of dnorm/2, 1-pval, (1-pval)/2

write.table(agg, file="TF_edges.tsv", sep="\t", quote=F, row.names=F)
