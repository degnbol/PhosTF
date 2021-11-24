#!/usr/bin/env zsh
# check if they are all PK
cut -f1 < scores.tsv | sed 1d | sort -u > unique_KPs.txt
grep $'\tPP' ../../network/KP_protein.tsv | cut -f1 | commul both - unique_KPs.txt
grep $'\tPK' ../../network/KP_protein.tsv | cut -f1 | commul both - unique_KPs.txt
# we see that they are all PK
