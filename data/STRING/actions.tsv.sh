#!/usr/bin/env zsh
# reversing protein1 and protein2 in interactions.tsv shows that all edges are given in both directions, so we can call protein1 source and protein2 target.
cat raw/4932.protein.actions.v*.txt | sed 's/4932\.//g' > actions.tsv
