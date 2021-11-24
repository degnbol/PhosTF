#!/usr/bin/env zsh
# some names related to tRNA has been started with capital T which is different from the standard lowercase. we correct.
sed $'/T.\(.*\)/s/\tT/\tt/' ../../raw/workman_2006/Suppl_Table_1.TF_Target.tsv > edges.tsv
# YGR040W aka. KSS1 is a protein kinase and it has edges in this dataset that was supposed to be for TF so it has been separated.
grep -v "^YGR040W" edges.tsv > TF_edges.tsv
grep "^YGR040W" edges.tsv > PK_edges.tsv
# add header
echo $'PK\tTarget\tSLL\tPval' | cat - PK_edges.tsv > temp && mv temp PK_edges.tsv
# average the duplicate entry: YNL309W YJL019W 
grep -v $'YNL309W\tYJL019W' TF_edges.tsv | cat - <(echo $'YNL309W\tYJL019W\t4.45223\t0.02451') > temp && mv temp TF_edges.tsv
