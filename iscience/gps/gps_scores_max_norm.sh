#!/usr/bin/env zsh
cat gps_scores_max.tsv | mlr --tsv stats1 -a sum -f score_max -g kinase then join -f gps_scores.tsv -j kinase then \
    put '$score_max_norm = $score_max / $score_max_sum' then cut -x -f score_max,score_max_sum | mlr --tsv sort -nr score_max_norm > gps_scores_max_norm.tsv
