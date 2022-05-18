#!/usr/bin/env zsh
./norm.R logFC_inner-strict.csv | sed '/1/s/""//' | gzip > logFC_inner-strict_norm.csv.gz
./norm.R logFC_outer-strict.csv | sed '/1/s/""//' | gzip > logFC_outer-strict_norm.csv.gz
