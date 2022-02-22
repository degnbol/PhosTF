#!/usr/bin/env zsh
../gene2ORF/gene2ORF.R edges_pop.tsv edges.tsv ../gene2ORF/gene2ORF.tsv
mlr --tsv --from edges.tsv filter '$KP2 == "" && $KP != $Target' then cut -x -f KP2 then uniq -a > KP_edges.tsv
