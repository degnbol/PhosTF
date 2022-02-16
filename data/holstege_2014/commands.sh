#!/usr/bin/env zsh
raw="../../raw/holstege_2014/deleteome_all_mutants_controls.txt"
cat $raw | cut -f2-3 | sed '2d' > gene2ORF.tsv
cat $raw | cut -f2 | sed '2d' | sed 's/systematicName/ORF/' > ORF.col
# number of columns
head $raw | awk -F'\t' '{print NF}' | sort | uniq
# 4464. there are 3*3 columns for wildtype at the end. A and M are swapped in order
cat $raw | cut -f2,4456- > wt.tsv
awk -F'\t' '{for(i=4;i<4456;i+=3){printf "\t%s", $i;} print "" }' $raw | sed '2d' > M.col
head -n1 M.col | tr '\t' '\n' | sed 's/-.*//' > header
# first convert using conversion builtin to rownames since we have previously tested that they are the same conversion results except for lug1 that is incorrectly referring to YCR087C-A instead of YLR352W in this data.
# I know to use YCR087C-A since the KO value is -3 instead of ~0 and since there would otherwise be two KOs of YLR352W, since it is present in the raw data.
../../../src/gene2ORF.R header header_ORF gene2ORF.tsv
../../../src/gene2ORF.R header_ORF header_ORF2 ../../name_conversion/gene2ORF.tsv
mv header_ORF2 header_ORF
# find lowercase entries, they have not been converted to anything
grep '[a-z]' header_ORF
# manually replace with the obvious names based on SGD searches
sed 's/yil014c/YIL014C-A/' header_ORF | sed 's/yol086w/YOL086W-A/' | sed 's/ydr034w/YDR034W-B/' | sed 's/yal044w/YAL044W-A/' > temp && mv temp header_ORF
cat header_ORF | tr '\n' '\t' | sed $'s/\t$//' > temp && mv temp header_ORF
sed '1d' M.col | cat header_ORF - | paste -d'\0' ORF.col - > M.tsv
rm header{,_ORF} *.col
# then remove duplicates with Rscript and maybe rm M.tsv afterwards to save space
