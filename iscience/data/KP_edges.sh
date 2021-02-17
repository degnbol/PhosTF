#!/usr/bin/env zsh
# a table of evaluation edges, and edges predicted with three methods: NetPhorest, GPS, and PhosTF. 
# It was filtered here for overlap but that reduced the data a lot and if we assume that any kinase that is modelled in the other methods are a random selection then we should not have to filter for comparison.
mlr --tsv --from ../../data/evaluation/P_eval.tsv rename Source,kinase,Target,substrate then filter '$netphorest != ""' then join -f ../gps/gps_edge_scores.tsv -j kinase,substrate then rename score_max_norm,gps then join -f ../../data/inference/KP_edges.tsv -j kinase,substrate > KP_edges.tsv
