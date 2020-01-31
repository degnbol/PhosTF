#!/usr/bin/env Rscript
library(Matrix)
library(ggplot2)

options(stringsAsFactors=FALSE)

# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/network")

# read
edges = read.table("TF_priors/TF_edges.tsv", sep="\t", header=T, quote="")
TFs = read.table("TF_mode.tsv", sep="\t", header=T, quote="")
Vs  = flatten(read.table("V.txt", quote=""))
simulated = flatten(read.table("../yeast_simulation/simulations/01/WT_est.mat", sep=" "))
simulated = simulated[simulated != 0]

# find a fit
plot(seq(0,1,.01), dbeta(seq(0,1,.01), 1, 20), type="l")
# create unsigned edge weights
alpha_ = 1
beta_ = 20
edges$weight = qbeta(edges$Pval, alpha_, beta_, lower.tail=F)

summary(edges$weight)
summary(abs(simulated))
sd(edges$weight)
# hist(edges$weight, breaks=200)
c(sum(edges$weight >= .1), sum(edges$weight >= .5), sum(edges$weight >= .75), sum(edges$weight >= .9))
mean(edges$weight[edges$weight > .75])

## add sign
edge_modes = TFs$Mode[match(edges$TF, TFs$TF)]
edges$Mode[edges$Mode == ""] = edge_modes[edges$Mode == ""]
# write repressor instead of both repressor and inhibitor
edges$Mode[edges$Mode == "inhibitor"] = "repressor"
# sum(edges$Mode == "activator")
# sum(edges$Mode == "repressor")
# edges$weight[edges$Mode == "activator"] = +edges$weight[edges$Mode == "activator"]  # redundant
edges$weight[edges$Mode == "repressor"] = -edges$weight[edges$Mode == "repressor"]

### plotting

plotdf = rbind(data.frame(weight=simulated, label="simulated"),
               data.frame(weight=edges$weight, label="constructed"))
labels = c("simulated", "constructed")
plotdf$label = factor(plotdf$label, levels=labels)
color_palette = c("gray50", "green4", "red")


bw = .05
plt = ggplot(plotdf, aes(weight, fill=label)) +
    geom_histogram(alpha=0.5, position="identity", binwidth=bw)
stepdf = ggplot_build(plt)$data[[1]][,c("xmin", "y", "group")]
stepdf$label = factor(stepdf$group, labels=labels)
# make visual corrections
for (i in 1:length(labels)) {
    limits = c(min(stepdf$xmin[stepdf$group==i])-bw, max(stepdf$xmin[stepdf$group==i])+bw)
    stepdf = rbind(stepdf, data.frame(xmin=limits, y=c(0,0), group=c(i,i), label=labels[i]))
}

pbreaks=c(-1e-60,-1e-12,-1e-6,-1e-3,-.05,.05,1e-3,1e-6,1e-12,1e-60)
px = sign(pbreaks) * qbeta(abs(pbreaks), alpha_, beta_, lower.tail=F)

plt + geom_step(data=stepdf, aes(x=xmin, y=y, color=label)) +
    theme_linedraw() +
    scale_color_manual(values=color_palette) +
    scale_fill_manual(values=color_palette) +
    theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
    ggtitle("TF edge weight construction") +
    ylab("edge count") +
    scale_y_log10(expand=c(0,0), limits=c(1,80000)) +
    scale_x_continuous(sec.axis=dup_axis(name="p-value", breaks=px, labels=abs(pbreaks)))


ggsave("edge_weights.pdf", width=7, height=2, units="in")


### write to files

write.table(edges, "TF_edge_weights.tsv", sep="\t", quote=F, row.names=F)

adjacency = as.matrix(sparseMatrix(i=match(edges$Target, Vs), j=match(edges$TF, TFs$TF), x=edges$weight, 
                                   dims=list(length(Vs), nrow(TFs)), dimnames=list(Vs, TFs$TF)))
adjacency_pval = as.matrix(sparseMatrix(i=match(edges$Target, Vs), j=match(edges$TF, TFs$TF), x=edges$Pval, 
                                   dims=list(length(Vs), nrow(TFs)), dimnames=list(Vs, TFs$TF)))
!any(edges$Pval==0)  # should be true then we can safely do:
adjacency_pval[adjacency_pval == 0] = NaN

adjacency_sign = sign(adjacency)
adjacency_sign[adjacency_sign == +1] = "+"
adjacency_sign[adjacency_sign == -1] = "-"

write.table(adjacency, "WT.csv", sep=",", quote=F)
write.table(adjacency, "WT.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacency_pval, "WT_p.csv", sep=",", quote=F, na="NaN")
write.table(adjacency_pval, "WT_p.mat", sep=" ", quote=F, col.names=F, row.names=F, na="NaN")
write.table(adjacency_sign, "WT_mask.csv", sep=",", quote=F)
write.table(adjacency_sign, "WT_mask.mat", sep=" ", quote=F, col.names=F, row.names=F)


noise_sd = 1/sqrt(sum(dim(adjacency)^2))
noise = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))
adjacency[adjacency == 0] = noise[adjacency == 0]

write.table(adjacency, "WT_noise.csv", sep=",", quote=F)
write.table(adjacency, "WT_noise.mat", sep=" ", quote=F, col.names=F, row.names=F)





