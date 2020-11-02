#!/usr/bin/env zsh
echo "norm is necessary since e.g. top 20% edges is all possible edges for some kinases while other kinases are simply filtered out fully."
cat gps_scores_max.tsv | mlr --tsv uniq -f kinase -c
head -n 40000 gps_scores_max.tsv | mlr --tsv uniq -f kinase -c
head -n 40000 gps_scores_max_norm.tsv | mlr --tsv uniq -f kinase -c
# we don't have a problem with substrates it seems
head -n 40000 gps_scores_max.tsv | mlr --tsv uniq -f substrate -c | wc -l
head -n 40000 gps_scores_max_norm.tsv | mlr --tsv uniq -f substrate -c | wc -l
