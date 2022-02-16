#!/usr/bin/env zsh
# running all against all
sed 1d ../biogrid/PTM.res > PTM.res
`git root`/etc/NetworKIN_v3/NetworKIN ../SGD/orfs.aa > networkin.tsv
`git root`/etc/NetworKIN_v3/NetworKIN ../SGD/orfs.aa PTM.res > networkin_biogrid.tsv
# reduce the heavy file by removing text description columns.
cat networkin.tsv | cut -f1-2,4-12,15-16 > temp && mv temp networkin.tsv
