#!/usr/bin/env zsh
# these are any undirected interactions in STRING
STRING="`git root`/data/STRING"
mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l protein1,protein2 -f $STRING/interactions.tsv then \
    filter '$EvalSet == 0' then cut -x -f EvalSet > fig4B_inSTRING.tsv
# has 7 entries

mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l item_id_a,item_id_b -f $STRING/actions.tsv > fig4B_actions.tsv
# has 5 entries

mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l Source,Target -f $STRING/directed.tsv > fig4B_directed.tsv
# has no entries

