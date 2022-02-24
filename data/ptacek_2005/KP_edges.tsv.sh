#!/usr/bin/env zsh
../gene2ORF/gene2ORF.R edges_pop.tsv edges.tsv ../gene2ORF/gene2ORF.tsv
# cyclin dependent kinases Pho85 and Cdc28 were measured with and without cyclins.
# Currently making the choice to discard cylin dependent edges.
mlr --tsv --from edges.tsv filter '$cyclin == "" && $KP != $Target' then cut -x -f KP2 then uniq -a > KP_edges.tsv
