#!/usr/bin/env zsh
mlr --tsv --from ../../data/evaluation/P_eval.tsv rename Source,kinase,Target,substrate then filter '$netphorest != ""' then join -f ../gps/gps_edge_scores.tsv -j kinase,substrate then rename score_max_norm,gps then join -f ../../data/inference/KP_edges.tsv -j kinase,substrate > KP_edges.tsv
