#!/usr/bin/env Rscript
library(Matrix)
library(ggplot2)

options(stringsAsFactors=FALSE)

# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/network")

# read
edges = read.table("TF_edges.tsv", sep="\t", header=T, quote="")
TFs = read.table("TF_mode.tsv", sep="\t", header=T, quote="")
Vs  = flatten(read.table("V.txt", quote=""))
simulated = flatten(read.table("../yeast_simulation/simulations/01/WT_est.mat", sep=" "))
simulated = abs(simulated[simulated != 0])

# create unsigned edge weights
edges$weight = qnorm(edges$Pval, mean=-.35, sd=.3, lower.tail=F)

plotdf = rbind(data.frame(weight=simulated[simulated>.1], label="simulated"),
               data.frame(weight=edges$weight, label="constructed"))
labels = c("simulated", "constructed")
plotdf$label = factor(plotdf$label, levels=labels)
color_palette = c("gray50", "green4", "red")

plt = ggplot(plotdf, aes(weight, y=..ncount.., fill=label)) +
    geom_histogram(alpha=0.5, position="identity", binwidth=.05)

stepdf = ggplot_build(plt)$data[[1]][,c("xmin", "y", "group")]
stepdf$label = factor(stepdf$group, labels=labels)
# make visual corrections
for (i in 1:length(labels)) {
    limits = c(min(stepdf$xmin[stepdf$group==i])-1e-2, max(stepdf$xmin[stepdf$group==i])+1e-2)
    stepdf = rbind(stepdf, data.frame(xmin=limits, y=c(0,0), group=c(i,i), label=labels[i]))
}

plt + geom_step(data=stepdf, aes(x=xmin, y=y, color=label)) +
    theme_linedraw() +
    scale_color_manual(values=color_palette) +
    scale_fill_manual(values=color_palette) +
    theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
    ggtitle("TF edge weights") +
    ylab("scaled count") +
    scale_x_continuous(sec.axis=sec_axis(~ pnorm(., mean=-.35, sd=.3, lower.tail=F), breaks=c(.05,1e-3,1e-6,1e-9,1e-15,1e-21), name="p-value"), limits=c(0,2.5))

ggsave("edge_weights.pdf", width=5, height=2.25, units="in")

summary(edges$weight)
sd(edges$weight)
# hist(edges$weight, breaks=200)
c(sum(edges$weight >= .5), sum(edges$weight >= .75), sum(edges$weight >= 1.))
mean(edges$weight[edges$weight > .75])


log_plot = function() {
    plt = ggplot(plotdf, aes(weight, fill=label)) +
        geom_histogram(alpha=0.5, position="identity", binwidth=.05)
    stepdf = ggplot_build(plt)$data[[1]][,c("xmin", "y", "group")]
    stepdf$label = factor(stepdf$group, labels=labels)
    # make visual corrections
    for (i in 1:length(labels)) {
        limits = c(min(stepdf$xmin[stepdf$group==i])-1e-2, max(stepdf$xmin[stepdf$group==i])+1e-2)
        stepdf = rbind(stepdf, data.frame(xmin=limits, y=c(0,0), group=c(i,i), label=labels[i]))
    }
    plt + geom_step(data=stepdf, aes(x=xmin, y=y, color=label)) +
        theme_linedraw() +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette) +
        theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
        ggtitle("TF edge weights") +
        ylab("scaled count") +
        scale_x_continuous(sec.axis=sec_axis(~ pnorm(., mean=-.35, sd=.3, lower.tail=F), breaks=c(.05,1e-3,1e-6,1e-9,1e-15,1e-21), name="p-value")) +
        scale_y_log10(expand=c(0,0))
}


edge_modes = TFs$Mode[match(edges$TF, TFs$TF)]
edges$Mode[edges$Mode == ""] = edge_modes[edges$Mode == ""]
# write repressor instead of both repressor and inhibitor
edges$Mode[edges$Mode == "inhibitor"] = "repressor"
# sum(edges$Mode == "activator")
# sum(edges$Mode == "repressor")

# edges$weight[edges$Mode == "activator"] = +edges$weight[edges$Mode == "activator"]  # redundant
edges$weight[edges$Mode == "repressor"] = -edges$weight[edges$Mode == "repressor"]

write.table(edges, "TF_edge_weights.tsv", sep="\t", quote=F, row.names=F)


adjacency = as.matrix(sparseMatrix(i=match(edges$Target, Vs), j=match(edges$TF, TFs$TF), x=edges$weight, 
                                   dims=list(length(Vs), nrow(TFs)), dimnames=list(Vs, TFs$TF)))
write.table(adjacency, "TF_edges.csv", sep=",", quote=F)
write.table(adjacency, "TF_edges.mat", sep=" ", quote=F, col.names=F, row.names=F)



