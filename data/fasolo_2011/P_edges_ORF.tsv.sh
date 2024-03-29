#!/usr/bin/env zsh
# correct systematic names by inserting a dash, e.g. "YHR214CC" -> "YHR214C-C"
sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])/\1-\2/' kinase_edges.tsv > kinase_edges-.tsv
../gene2ORF/gene2ORF.R kinase_edges-.tsv P_edges_ORF.tsv ../gene2ORF/gene2ORF.tsv
rm kinase_edges-.tsv

