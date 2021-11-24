awk '{print $2 "\t" $1}' ../../raw/chua_2006/{Del,OE}_NormalizedRatios.phenoData | grep -v ORF | sed $'/-$/s/\t/-\t/' | sed $'/+$/s/\t/+\t/' > gene2ORF.tsv
head -n1 ../../raw/chua_2006/OE_NormalizedRatios.txt | tr -d '\r' > OE.head
head -n1 ../../raw/chua_2006/Del_NormalizedRatios.txt | tr -d '\r' > KO.head
# conversion of names was checked to be the same as if it was done with the converter in name_conversion folder.
../../../src/gene2ORF.R OE.head OE_ORF.head gene2ORF.tsv
r OE=KO
# exp column of .PhenoData was found to be in the same order as the header of the measurement file, so the ORF column could be used for replacing the header directly but I used conversion tool anyways.
# micro array naming error with naming things ending in 01 and the dash is removed. Corrected:
cat ../../raw/chua_2006/OE_NormalizedRatios.txt | tr -d '\r' | sed 1d | sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])01/\1-\2/' | cat OE_ORF.head - > OE.tsv
cat ../../raw/chua_2006/Del_NormalizedRatios.txt | tr -d '\r' | sed 1d | sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])01/\1-\2/' | cat KO_ORF.head - > KO.tsv
# check that signs do not need to be reversed on any measurements
Rscript check_signs.R
# remove sign notations
cat OE_ORF.head | tr -d '-' | tr -d '+' > OE_ORF_uns.head
r OE=KO
sed 1d OE.tsv | cat OE_ORF_uns.head - > TF_OE.tsv
r OE=KO
rm *.head
