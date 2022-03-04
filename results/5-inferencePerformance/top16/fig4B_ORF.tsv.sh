#!/usr/bin/env zsh
# REQUIRES: fig4B.tsv, gene2ORF.tsv gene2ORF.R
G2O="`git root`/data/gene2ORF"
$G2O/gene2ORF.R fig4B.tsv fig4B_ORF.tsv $G2O/gene2ORF.tsv
mlr --tsv --from fig4B.tsv rename P,P_gene,Target,Target_gene then \
    cut -x -f EvalSet | paste - fig4B_ORF.tsv > temp && mv temp fig4B_ORF.tsv
