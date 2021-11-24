#!/usr/bin/env zsh
head -n1 PK_KO.tsv | tr '\t' '\n' | sed 1d | sort -u > KPs.txt
