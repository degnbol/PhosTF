# there is some text at the bottom after some empty lines so we cut off at appropriate length with head
sed 1d ../../raw/harbison_2004/pvalbygene_forpaper_abbr_YPD.txt | tr -d '\r' | head -n6290 | sed -E '/^[[:space:]]+$/d' | sed 's/_YPD//g' > YPD.tmp
sed 1d ../../raw/harbison_2004/pvalbygene_forpaper_abbr_other_conditions.txt | tr -d '\r' | sed -E '/^[[:space:]]+$/d' | cut -d$'\t' -f1,4- > conds.tsv
sed 1d YPD.tmp | cut -d$'\t' -f1,2 > gene2ORF.tsv
cat YPD.tmp | cut -d$'\t' -f1,4- > YPD.tsv
rm YPD.tmp
# it was tested that the built in conversion matched perfectly in name conversion with the general one, but has less matched for the standard names to be converted.
head -n1 YPD.tsv | awk '{print toupper($0)}' > YPD.head
head -n1 conds.tsv | sed -E 's/_[^[:space:]]+//g' > conds.head
../../../src/gene2ORF.R YPD.head YPD_ORF.head ../../name_conversion/gene2ORF.tsv
r YPD=conds
sed 's/A1 (MATA1)/MATA1/' YPD_ORF.head | cat - <(sed 1d YPD.tsv) > temp && mv temp YPD.tsv
sed 1d conds.tsv | cat conds_ORF.head - > temp && mv temp conds.tsv
rm *.head 
# NOTE that there are 3 kinases that has been analysed here as if they were TFs, e.g. KSS1.
