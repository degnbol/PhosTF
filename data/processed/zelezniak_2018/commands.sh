#!/usr/bin/env zsh
echo -n $'ORF\t' | cat - PK_KO.tsv > temp && mv temp PK_KO.tsv
head -n1 PK_KO.tsv | tr '\t' '\n' | sed 1d | sort -u > KPs.txt
