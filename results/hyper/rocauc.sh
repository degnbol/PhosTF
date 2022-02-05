#!/usr/bin/env zsh
# $1 should be $hyper, $2 should be $λBstar-$λabsW
hyper=$1
reg=$2
outfile="aucs/auc-${1}-${2}.txt"
touch $outfile
for n in 10 100; do
    for i in {1..5}; do
        for rep in {1..5}; do
            echo -n "$n\t$i\t$rep\t$hyper\t$reg\tWT\t" >> $outfile
            ./rocauc.py "../2-insilicoNetworkSimulation/adjacencies/WT_${n}_${i}-rep${rep}.adj" "inferredWeights/$1/WT_${n}_${i}-rep${rep}-$2.adj" | grep -v '^#' >> $outfile
            echo -n "$n\t$i\t$rep\t$hyper\t$reg\tWP\t" >> $outfile
            ./rocauc.py "../2-insilicoNetworkSimulation/adjacencies/WP_${n}_${i}-rep${rep}.adj" "inferredWeights/$1/WP_${n}_${i}-rep${rep}-$2.adj" | grep -v '^#' >> $outfile
        done
    done
done
