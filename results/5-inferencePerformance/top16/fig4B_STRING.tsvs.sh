#!/usr/bin/env zsh
# REQUIRES: fig4B_ORF.tsv, STRING files interactions.tsv, actions.tsv, directed.tsv
# WRITES: fig4B_{inSTRING,inSTRING-score500,actions,directed,undirected-ptmod}.tsv
STRING="`git root`/data/STRING"
# since undirected we match against protein1 and protein2 as either source or target
mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l protein1,protein2 -f $STRING/interactions.tsv then \
    filter '$EvalSet == 0' then cut -x -f EvalSet > fig4B_inSTRING.tsv
mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l protein2,protein1 -f $STRING/interactions.tsv then \
    filter '$EvalSet == 0' then cut -x -f EvalSet | sed 1d >> fig4B_inSTRING.tsv
# has 14 entries
mlr --tsv --from fig4B_inSTRING.tsv filter '$combined_score > 500' > fig4B_inSTRING-score500.tsv
# has 10 entries

# since undirected we match against protein1 and protein2 as either source or target
mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l item_id_a,item_id_b -f $STRING/actions.tsv > fig4B_actions.tsv
# mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l item_id_b,item_id_a -f $STRING/actions.tsv | sed 1d >> fig4B_actions.tsv
# has 5 non-unique entries. Switching item_id_a and b gives duplicates so both directions are alreadt listed in $STRING/actions.tsv

mlr --tsv --from fig4B_ORF.tsv join -j P,target -r P,Target -l Source,Target -f $STRING/directed.tsv > fig4B_directed.tsv
# has no entries

# however I realize there are two undirected ptmod entry matches, one from eval set and one novel: 
mlr --tsv --from fig4B_actions.tsv filter '$mode == "ptmod"' > fig4B_undirected_ptmod.tsv


