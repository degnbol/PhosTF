#!/usr/bin/env zsh
# combines experimental data describing evidence of kinase and phosphatase target relationships.
# REQUIRES: multiple processed data files
DATA=`git root`/data

mlr --tsv --from $DATA/STRING/scores.tsv filter '$ptmod != ""' then \
    uniq -f Source,Target,ptmod then rename Source,P,ptmod,Value then \
    put '$Variable = "ptmod"; $Ref = "STRING"' > P_data-ptmod.tsv

mlr --tsv --from $DATA/biogrid/P_edges.tsv rename Relationship,Value then \
    put '$Variable = "Relationship"; $Ref = "BioGRID"' > P_data-biogrid.tsv

# no value associated with each edge, only that they are present in the list of edges.
mlr --tsv --from $DATA/fasolo_2011/P_edges_ORF.tsv rename kinase,P,target,Target then \
    put '$Value = ""; $Variable = ""; $Ref = "Fasolo 2011"' > P_data-fasolo.tsv

mlr --tsv --from $DATA/yeastkid/yeastkid.tsv rename 'Kinase ID,P,Gene ID,Target,Score,Value' then \
    uniq -f P,Target,Value then \
    put '$Variable = "Score"; $Ref = "Yeast KID"' > P_data-yeastkid.tsv

# no value associated with each edge, only that they are present in the list of edges.
mlr --tsv --from $DATA/ptacek_2005/KP_edges.tsv rename KP,P then \
    put '$Value = ""; $Variable = ""; $Ref = "Ptacek 2005"' > P_data-ptacek.tsv

mlr --tsv cat P_data-*.tsv > P_data.tsv

