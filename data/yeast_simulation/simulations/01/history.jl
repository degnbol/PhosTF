#!/usr/bin/env julia
# in julia interactive. _cor files made with
PKTFX.correct("WT_rand.mat")
PKTFX.network("WT_cor.mat", "WP_cor.mat")
PKTFX.estimateWt(; o="WT_est.mat")
