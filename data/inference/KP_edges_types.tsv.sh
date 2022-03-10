#!/usr/bin/env zsh
mlr --tsv --from ../network/node_attributes_full.txt cut -f ORF,protein_type then\
    join -j Target -l Target -r ORF -f KP_edges.tsv then rename protein_type,Target_type then\
    join -j P -r KP -l ORF -f ../network/KP_protein.tsv then rename protein_type,P_type then\
    reorder -f P,P_type,Target,Target_type then\
    sort -n q > KP_edges_types.tsv

# annotation files ../network/node_attributes_full.txt and ../network/KP_protein.tsv should have unique entries
# Otherwise we would be doing a cartesian join.
if [ `wc -l < KP_edges.tsv` -ne `wc -l < KP_edges_types.tsv` ]; then
    echo "Join problem"
    rm KP_edges_types.tsv
fi

