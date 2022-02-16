#!/usr/bin/env zsh
# All non yeast entries removed and redundant columns as well.
# sort | uniq performed.
grep Phosphorylation raw/BIOGRID-PTM-3.5.178.ptmtab.tsv | cut -f4 | sed 1d | sort -u > targets.txt
cat raw/BIOGRID-PTM-3.5.178.ptmtab.tsv | cut -f4,9,11 | grep -v $'\t-' | sort -u > PTM.res
cat raw/BIOGRID-PTM-3.5.178.ptmtab.tsv | cut -f4,7 | sed 1d | grep -v $'\t-' | sort -u > sequences.linfa
cat sequences.linfa | sed 's/^/>/' | tr '\t' '\n' > sequences.fa
grep kinase raw/BIOGRID-PTM-RELATIONSHIPS-3.5.178.ptmrel.tsv | cut -f4 | sort -u > kinases.txt
r kinase=phosphatase
cat <(sed $'s/$/\tkinase/' kinases.txt) <(sed $'s/$/\tphosphatase/' phosphatases.txt) > KP.tsv
# number of unique sites
cat raw/BIOGRID-PTM-RELATIONSHIPS-3.5.178.ptmrel.tsv | cut -f1,7 | grep -v $'\t-' | sed 1d | cut -f1 | sort -n | uniq | wc -l
