#!/usr/bin/env zsh
gzcat raw/table.txt.gz | sed '1,8d' | mlr --tsvlite uniq -f 'Gene Name,Kinase Name,Gene ID,Kinase ID,Score' | gzip > yeastkid.tsv.gz
