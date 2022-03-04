#!/usr/bin/env zsh
# all 12 have ptmod values given
mlr --tsv --from fig4B_inSTRING.tsv filter '$ptmod != ""' then stats1 -a max,mean -f ptmod -g P_gene,Target_gene > fig4B_ptmod.tsv
