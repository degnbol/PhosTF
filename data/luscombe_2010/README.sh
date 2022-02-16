#!/usr/bin/env zsh
# adding ORF at beginning of file to be consistent with style from other filesets.
echo -n "ORF\t" | cat - TF_KO.tsv > temp && mv temp TF_KO.tsv
r TF=PK
r KO=KO_pval
r PK=TF
