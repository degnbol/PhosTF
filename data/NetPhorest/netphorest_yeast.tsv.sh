#!/usr/bin/env zsh
cat raw/Saccharomyces_cerevisiae.EF4.74.pep.all.fa | `git root`/etc/NetPhorest_yeast_v2/netphorest > netphorest.out
# an extra column name was inserted for KP since it was apparently missing
grep yeast netphorest.tsv | cat <(head -n1 netphorest.tsv | sed 's/# //') - > netphorest_yeast.tsv
cat netphorest_yeast.tsv | cut -f1,8,9 | grep -v group > netphorest_yeast_col189.tsv
