#!/usr/bin/env zsh
./commands.jl
# for mutnum in {1..4}; do
# `git root`/src/simulation/plot_simulated_mutant.R 2 2 1 $mutnum
# done
`git root`/src/inference/infer.jl sim_logFC.tsv 'T.*' 'P.*' --epochs 30000
./xgmmls.jl
