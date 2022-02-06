#!/usr/bin/env zsh
./print_best.sh | head -n2 | mlr --tsv join -f aucs.tsv -j cancel,mean_k,vec,norm,Bstar,absW then sort -n AUC
