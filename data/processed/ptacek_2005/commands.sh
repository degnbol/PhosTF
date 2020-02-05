#!/usr/bin/env zsh
../../../src/gene2ORF.R edges_pop.tsv edges.tsv ../../name_conversion/gene2ORF.tsv
grep '\t\t' edges.tsv | sed $'s/\t\t/\t/' | cat <(echo $'KP\tTarget') - > KP_edges.tsv
