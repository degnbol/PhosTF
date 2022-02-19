#!/usr/bin/env zsh
for iORo in "inner" "outer"; do
    for logFC in "-strict" "-enh1" "-enh4"; do
        for lB in "0.0" "0.01" "0.1" "1.0"; do
            for lW in "0.0" "0.01" "0.1" "1.0"; do
                ./infer.jl $iORo $logFC $lB $lW > log$iORo$logFC-B$lB-W$lW.tsv &
            done
        done
    done
done
