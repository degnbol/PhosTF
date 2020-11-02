#!/usr/bin/env zsh
cat gps_scores.tsv | mlr --tsv stats1 -a sum -f score -g kinase then join -f gps_scores.tsv -j kinase then \
    put '$score_norm = $score / $score_sum' then cut -x -f score,score_sum | mlr --tsv sort -nr score_norm > gps_scores_norm.tsv
