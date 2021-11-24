#!/usr/bin/env zsh
echo -n $'ORF\t' | cat - PK_KO.tsv > temp && mv temp PK_KO.tsv
