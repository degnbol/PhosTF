#!/usr/bin/env zsh
cat ../../raw/NetPhorest/Saccharomyces_cerevisiae.EF4.74.pep.all.fa | ../../../apps/NetPhorest_yeast_v2/netphorest > netphorest.out
# an extra column name was inserted for KP since it was apparently missing
grep yeast netphorest.tsv | cat <(head -n1 netphorest.tsv | sed 's/# //') - > netphorest_yeast.tsv
cat netphorest_yeast.tsv | cut -d$'\t' -f1,8,9 | grep -v group > netphorest_yeast_col189.tsv
# check if they are all PK
cut -f1 < scores.tsv | sed 1d | sort -u > unique_KPs.txt
grep $'\tPP' ../../network/KP_protein.tsv | cut -f1 | commul both - unique_KPs.txt
grep $'\tPK' ../../network/KP_protein.tsv | cut -f1 | commul both - unique_KPs.txt
# we see that they are all PK
