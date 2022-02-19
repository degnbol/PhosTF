#!/usr/bin/env zsh
for iORo in "inner" "outer"; do
    for logFC in "" "-strict" "-enh1" "-enh4"; do
        for lB in "0.0" "0.01" "0.1" "0.5" "1.0"; do
            for lW in "0.0" "0.01" "0.1" "0.5" "1.0"; do
                julia -t 10 ./infer.jl $iORo $logFC $lB $lW &
            done
        done
    done
done
