#!/usr/bin/env zsh
cut -f 1,3,5 yeast.gps | sed $'1s/$/\tsubstrate/' | sed $'/^>/s/>/>\t\t\t/' | awk -F$'\t' 'BEGIN{OFS="\t"} $4=="" {$4=prev} {prev=$4}1' | grep -v '^>' |
    cut -f2- | mlr --tsv stats1 -a max -f Score -g Kinase,substrate then rename Kinase,kinase,Score_max,score > gps_scores.tsv
