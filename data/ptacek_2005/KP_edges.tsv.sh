#!/usr/bin/env zsh
`git root`/src/gene2ORF.R edges_pop.tsv edges.tsv `git root`/src/gene2ORF.tsv
grep '\t\t' edges.tsv | sed $'s/\t\t/\t/' | cat <(echo "KP\tTarget") - | awk '$1 != $2' | sort | uniq > KP_edges.tsv
