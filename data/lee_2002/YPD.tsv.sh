#!/usr/bin/env zsh
# pvals for TF gene binding
sed 1d raw/lee_2002/106_pvalbygene_ypd_v9.0.txt | tr -d '\r' | head -n6290 | sed -E '/^[[:space:]]+$/d' | cut -f1,4- > YPD.tsv
head -n1 YPD.tsv > YPD.head
../gene2ORF/gene2ORF.R YPD.head YPD_ORF.head ../gene2ORF/gene2ORF.tsv
sed 's/A1 (MATA1)/MATA1/' YPD_ORF.head | cat - <(sed 1d YPD.tsv) > temp && mv temp YPD.tsv
rm *.head
