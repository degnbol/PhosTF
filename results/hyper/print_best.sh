#!/usr/bin/env zsh
mlr --tsv --from aucs.tsv stats1 -a mean -f AUC -g cancel,mean_k,vec,norm,Bstar,absW then sort -nr AUC_mean | head
