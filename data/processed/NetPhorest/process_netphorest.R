
setwd("~/cwd/data/processed/NetPhorest")

posteriors = read.table("netphorest_yeast_col189.tsv", sep="\t", header=T, quote="")
posteriors_agg = aggregate(Posterior ~ ., data=posteriors, max)
scores = data.frame(KP=posteriors_agg$KP, Target=posteriors_agg$Name, Score=posteriors_agg$Posterior)
write.table(scores, file="scores.tsv", sep="\t", quote=F, row.names=F)
# scores = read.table("scores.tsv", sep="\t", quote="", header=T)

plot(density(scores$Score))
hist(scores$Score, breaks=30)
sum(scores$Score > 0.3)
sum(scores$Score < 0.03)
