# pvals for TF gene binding
sed 1d ../../raw/lee_2002/106_pvalbygene_ypd_v9.0.txt | tr -d '\r' | head -n6290 | sed -E '/^[[:space:]]+$/d' | cut -d$'\t' -f1,4- > YPD.tsv
head -n1 YPD.tsv > YPD.head
../../../src/gene2ORF.R YPD.head YPD_ORF.head ../../name_conversion/gene2ORF.tsv
sed 's/A1 (MATA1)/MATA1/' YPD_ORF.head | cat - <(sed 1d YPD.tsv) > temp && mv temp YPD.tsv
rm *.head
# lee 2002 is from same lab (young lab http://younglab.wi.mit.edu/) as harbison 2004 and the data is the same for a lot of the entries.
# we combine both in harbison, so this data is used there and should not be used on its own.
