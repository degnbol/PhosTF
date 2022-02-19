#!/usr/bin/env Rscript

library(ggplot2)

# function
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/processed/holstege_2010")
# read
PKKO = read.table("PK_KO_Mpval.tsv", header=T, sep="\t", quote="", stringsAsFactors=F)

M = abs(flatten(PKKO[,seq(3,ncol(PKKO),2)]))
p = flatten(PKKO[,seq(4,ncol(PKKO),2)])

df = data.frame(M=M[p < 0.05], p=p[p < 0.05])
df = df[df$p != 0,]
df$logp = log10(df$p)
df$logM = log2(df$M)

ggplot(df, aes(M,p)) +
    geom_point(alpha=0.01) +
    scale_x_continuous(trans="log2") +
    scale_y_log10()
 
ggplot(df, aes(logM,logp)) +
    geom_point(alpha=0.01)

# lm(logp ~ logM, df)
fit = lm(y ~ x, data.frame(x=c(-2.5, -1.5), y=c(-1,-6)))
coefs = fit$coefficients

# it fits
mean(df$logp - (coefs[2] * df$logM + coefs[1]))
M2p = function(M) as.vector(10^(coefs[2] * log2(M) + coefs[1]))
p2M = function(p) as.vector(2^((log10(p) - coefs[1]) / coefs[2]))

p2M(c(0.05,0.001))
M2p(1:5)

df$fitp = M2p(df$M)
ggplot(df) +
    geom_point(aes(M,p), alpha=0.01) +
    geom_line(aes(M,fitp), color="red") +
    scale_x_continuous(trans="log2") +
    scale_y_log10(limits=c(1e-07,1e-01))

ggplot(df) +
    geom_point(aes(M,p), alpha=0.01) +
    geom_line(aes(M,fitp), color="red") +
    ylim(0,0.05)








