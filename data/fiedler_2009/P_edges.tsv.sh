#!/usr/bin/env zsh
# remove ending commas so we can select all rows with 5 fields. Remove two headers with grep.
sed 's/,*$//' raw/mmc4.csv | awk -F, "NF==5" | mlr --icsv --otsv cat |
    grep -v 'Type of interaction' | 
    cat <(echo "P_ORF\tP_gene\ttarget_ORF\ttarget_gene\tmode") - > P_edges.tsv
# NOTE: ignores the last part of the file where kinase and phosphatase dependence is described.
