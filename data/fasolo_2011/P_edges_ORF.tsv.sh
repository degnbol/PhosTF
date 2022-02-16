#!/usr/bin/env zsh
# correct sys names by inserting a -
sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])/\1-\2/' P_edges_ORF.tsv > temp && mv temp P_edges_ORF.tsv
