#!/usr/bin/env zsh
mlr --tsv --from aucs.tsv filter '$source == "P"' then stats1 -a mean,min -f AUC -g cancel,mean_k,vec,norm,Bstar,absW then sort -nr AUC_mean,AUC_min | head
