#!/usr/bin/env zsh
mlr --icsv --otsv --from raw/mmc2.csv cut -x -f 'SGD annotation' then uniq -a > geneType.tsv
