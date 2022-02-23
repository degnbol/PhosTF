#!/usr/bin/env zsh
# REQUIRES: P_data-biogrid.tsv and P_edges.tsv
# INSTALL: merge.R from github.com/degnbol/degnlib
# annotate P_edges.tsv with the kinase or phosphatase action from biogrid
merge.R -H P_data-biogrid.tsv < P_edges.tsv | 
    mlr --tsv cut -f P,Target,Value then rename Value,Mode > P_edges_mode.tsv
