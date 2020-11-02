#!/usr/bin/env zsh
# All tyrosine kinases seem to be under the category Dual but at the same time the next groupings (separated with /) does not seem to be a protein (one is e.g. called "other")
# We try to filter out values here that are for kinase families and only keep values for single kinases.
mlr --tsv --from gps_scores.tsv put '$kinase = sub($kinase, "Dual/", "")' then filter '$kinase =~ "/"' then put '$kinase = sub($kinase, ".*/", "")' then \
    stats1 -a max -f score -g kinase,substrate | mlr --tsv sort -nr score_max > gps_scores_max.tsv
