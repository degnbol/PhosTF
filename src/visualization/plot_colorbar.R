#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr) # so we can take a colorbar from a plot and only plot that
suppressPackageStartupMessages(library(here))

# load precompiled color map
rgbs = read.table(paste0(here(), "/src/utilities/colormap_D7.txt"))
len = ncol(rgbs)
alpha = abs(4 * ((1:len -1) / len - .5)); alpha[alpha >= 1] = 1
cmap = rgb(rgbs[1,], rgbs[2,], rgbs[3,])
cmap_alpha = rgb(rgbs[1,], rgbs[2,], rgbs[3,], alpha)

df = data.frame(x=seq(-1, 1, length.out=len), y=0)

p0 = ggplot(df, aes(x=x, y=y, color=x)) +
    geom_line(size=10) +
    xlim(-1, 1) + ylim(0,001) +
    theme_minimal() + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

guid = guide_colourbar(title="Node value", barwidth=1, barheight=10, title.position="bottom", nbin=256)
p = p0 + xlab("Node value") + scale_colour_gradientn(colors=cmap, guide=guid, breaks=seq(-1,1,.25))
as_ggplot(get_legend(p))

guid = guide_colourbar(title="Edge value", barwidth=1, barheight=10, title.position="bottom", nbin=256)
p = p0 + xlab("Edge value") + scale_colour_gradientn(colors=cmap_alpha, guide=guid, breaks=seq(-1,1,.25))
as_ggplot(get_legend(p))

