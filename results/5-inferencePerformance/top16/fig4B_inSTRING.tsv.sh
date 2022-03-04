#!/usr/bin/env zsh
mlr --tsv --from fig4B_ORF.tsv filter '$EvalSet == 0' then cut -x -f EvalSet then join -j P -l Source -r P -f `git root`/data/STRING/scores.tsv > fig4B_inSTRING.tsv
