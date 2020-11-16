#!/usr/bin/env zsh
mlr --tsv rename KP,kinase,Target,substrate 03/KP_edges.tsv | mlr --tsv sort -n q > KP_edges.tsv
