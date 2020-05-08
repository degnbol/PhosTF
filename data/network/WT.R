#!/usr/bin/env Rscript
library(Matrix)
library(ggplot2)
library(fdrtool)

options(stringsAsFactors=FALSE)

# functions
flatten = function(x) as.vector(as.matrix(x))
qhalfnorm_ = function(pvals, sd) {
    vals = qhalfnorm(edges$Pval, theta=sd2theta(sd), lower.tail=F)
    vals[vals == Inf] = max(vals[vals != Inf])
    vals
}
sparsematrix = function(i, j, x) {
    i = match(i, Vs)
    j = match(j, TFs$TF)
    dims_ = list(length(Vs), nrow(TFs))
    dimnames_ = list(Vs, TFs$TF)
    as.matrix(sparseMatrix(i=i, j=j, x=x, dims=dims_, dimnames=dimnames_))
}

add_noise = function(adjacency, noise_sd = 1/sqrt(prod(dim(adjacency)))) {
    cat(noise_sd, "\n")
    adjacency_noise = adjacency
    lacking = adjacency_noise==0
    adjacency_noise[lacking] = matrix(rnorm(prod(dim(adjacency)), sd=noise_sd), nrow=nrow(adjacency), ncol=ncol(adjacency))[lacking]
    adjacency_noise
}

get_sign_mask = function(adjacency) {
    out = sign(adjacency)
    out[out == +1] = "+"
    out[out == -1] = "-"
    out
}


setwd("~/cwd/data/network")

# read
edges = read.table("TF_priors/TF_edges.tsv", sep="\t", header=T, quote="")
TFs = read.table("TF_mode.tsv", sep="\t", header=T, quote="")
Vs  = flatten(read.table("V.txt", quote=""))
simulated = flatten(read.table("../yeast_simulation/simulations/01/WT_est.mat", sep=" "))
simulated = simulated[simulated != 0]
nTF = nrow(TFs)
nV = length(Vs)

# find a fit
plot(seq(0,1,.01), dbeta(seq(0,1,.01), 1, 20), type="l")
# create unsigned edge weights
edges$gauss01 = qhalfnorm_(edges$Pval, 0.10)
edges$gauss025 = qhalfnorm_(edges$Pval, 0.25)
edges$weight = qbeta(edges$qval, 1, 20, lower.tail=F)
edges$weight25 = qbeta(edges$qval, 1, 25, lower.tail=F)


summary(edges$weight25)
summary(abs(simulated))
sd(edges$weight25)
hist(edges$weight25, breaks=200)
c(sum(edges$weight25 >= .1), sum(edges$weight25 >= .5), sum(edges$weight25 >= .75), sum(edges$weight25 >= .9))
mean(edges$weight25[edges$weight25 > .75])

## add sign
edge_modes = TFs$Mode[match(edges$TF, TFs$TF)]
edges$Mode[edges$Mode == ""] = edge_modes[edges$Mode == ""]
# write repressor instead of both repressor and inhibitor
edges$Mode[edges$Mode == "inhibitor"] = "repressor"
# sum(edges$Mode == "activator")
# sum(edges$Mode == "repressor")
edges$weight[edges$Mode == "repressor"] = -edges$weight[edges$Mode == "repressor"]
edges$weight25[edges$Mode == "repressor"] = -edges$weight25[edges$Mode == "repressor"]

edges_FDR10 = edges[edges$qval<.1,]
edges_FDR20 = edges[edges$qval<.2,]

### plotting

weight_plot = function(dataset, beta_) {
    plotdf = rbind(data.frame(weight=simulated, label="simulated"),
                   data.frame(weight=dataset, label="constructed"))
    labels = c("simulated", "constructed")
    plotdf$label = factor(plotdf$label, levels=labels)
    color_palette = c("gray50", "green4", "red")
    bw = .05
    
    Beta_curve = data.frame(weight=seq(-1,1,bw/2), dens=dbeta(abs(seq(-1,1,bw/2)), 1, beta_))
    
    plt = ggplot(plotdf, aes(weight)) +
        geom_histogram(mapping=aes(fill=label), alpha=0.5, position="identity", binwidth=bw) #+
        #geom_line(data=Beta_curve, mapping=aes(y=2500*dens+1), linetype=2)
    
    
    stepdf = ggplot_build(plt)$data[[1]][,c("xmin", "y", "group")]
    stepdf$label = factor(stepdf$group, labels=labels)
    # make visual corrections
    for (i in 1:length(labels)) {
        limits = c(-1-bw,1+bw/2)
        stepdf = rbind(stepdf, data.frame(xmin=limits, y=c(0,0), group=c(i,i), label=labels[i]))
    }
    
    pbreaks=c(-1e-60,-1e-12,-1e-6,-1e-3,-.05,.05,1e-3,1e-6,1e-12,1e-60)
    px = sign(pbreaks) * qbeta(abs(pbreaks), 1, beta_, lower.tail=F)
    
    plt + geom_step(data=stepdf, aes(x=xmin, y=y, color=label)) +
        theme_linedraw() +
        scale_color_manual(values=color_palette) +
        scale_fill_manual(values=color_palette) +
        theme(legend.title=element_blank(), panel.grid.major=element_line(colour="lightgray"), panel.grid.minor=element_blank()) +
        ggtitle("TF edge weight construction") +
        ylab("edge count") +
        scale_y_log10(expand=c(0,0), limits=c(1,80000)) +
        scale_x_continuous(sec.axis=dup_axis(name="p-value", breaks=px, labels=abs(pbreaks)), limits=c(-1.1-bw,1.1+bw), expand=c(0,0))
        
}


weight_plot(edges$weight, 20)
weight_plot(edges_FDR10$weight, 20)
weight_plot(edges_FDR20$weight, 20)
# ggsave("edge_weights.pdf", plot=weight_plot(edges_FDR10$weight25, 25), width=7, height=2, units="in")


### write to files

write.table(edges_FDR20, "TF_edge_weights.tsv", sep="\t", quote=F, row.names=F)

get_adjacencies = function(edges, weight) {
    adjacency = sparsematrix(edges$Target, edges$TF, edges$weight)
    adjacency_pval = sparsematrix(edges$Target, edges$TF, edges$Pval)
    stopifnot(!any(edges$Pval==0))  # should be true then we can safely do:
    adjacency_pval[adjacency_pval == 0] = NaN
    list(adjacency, adjacency_pval, get_sign_mask(adjacency), add_noise(adjacency))
}

adjacencies = get_adjacencies(edges, "weight")
write.table(adjacencies[[1]], "WT.csv", sep=",", quote=F)
write.table(adjacencies[[1]], "WT.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacencies[[2]], "WT_p.csv", sep=",", quote=F, na="NaN")
write.table(adjacencies[[2]], "WT_p.mat", sep=" ", quote=F, col.names=F, row.names=F, na="NaN")
write.table(adjacencies[[3]], "WT_mask.csv", sep=",", quote=F)
write.table(adjacencies[[3]], "WT_mask.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(adjacencies[[4]], "WT_noise.csv", sep=",", quote=F)
write.table(adjacencies[[4]], "WT_noise.mat", sep=" ", quote=F, col.names=F, row.names=F)

adjacencies_FDR10 = get_adjacencies(edges_FDR10, "weight25")
write.table(adjacencies_FDR10[[1]], "WT_FDR10.csv", sep=",", quote=F)
write.table(adjacencies_FDR10[[1]], "WT_FDR10.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(sign(adjacencies_FDR10[[1]]), "WT_FDR10_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)

adjacencies_FDR20 = get_adjacencies(edges_FDR20, "weight25")
write.table(adjacencies_FDR20[[1]], "WT_FDR20.csv", sep=",", quote=F)
write.table(adjacencies_FDR20[[1]], "WT_FDR20.mat", sep=" ", quote=F, col.names=F, row.names=F)
write.table(sign(adjacencies_FDR20[[1]]), "WT_FDR20_sign.mat", sep=" ", quote=F, col.names=F, row.names=F)


fwrite(sparsematrix(edges$Target, edges$TF, edges$gauss01), "WT_gauss01.mat", sep=" ", row.names=F, col.names=F)








