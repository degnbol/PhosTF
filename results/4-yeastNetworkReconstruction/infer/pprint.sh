#!/usr/bin/env zsh
tr '+' '\t' < hyperps.tsv |
    sed 's/.tsv//' | sed 's/WP_infer-//' | sed 's/ $//' |
    cat <(echo "DF_size\tDF_issues\tB*\t|W|\tWT0\tWT_mask\tp") - |
    mlr --itsv --opprint sort -n p
