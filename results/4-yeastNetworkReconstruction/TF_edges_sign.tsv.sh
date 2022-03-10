#!/usr/bin/env zsh
`git root`/src/scripts/melt.R WT_FDR20_sign_mask.ssv Target TF Sign |
    mlr --tsv filter '$Sign != "."' then reorder -f TF > TF_edges_sign.tsv
