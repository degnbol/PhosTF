#!/usr/bin/env zsh
# annotate P_edges.tsv with the kinase or phosphatase action from biogrid
merge.R -H P_data-biogrid.tsv < P_edges.tsv | 
    mlr --tsv cut -f P,Target,Value then rename Value,Mode > P_edges_mode.tsv
