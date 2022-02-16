#!/usr/bin/env zsh
# prepend ORF
echo -n "ORF\t" | cat - PK_KO.tsv > temp && mv temp PK_KO.tsv
