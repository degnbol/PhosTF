#!/usr/bin/env zsh
mlr --tsv --from P_eval.tsv cut -f Source,Target,netphorest then filter '$netphorest != ""' | mlr --tsv sort -nr netphorest > netphorest.tsv
