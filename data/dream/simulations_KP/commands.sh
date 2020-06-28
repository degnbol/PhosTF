#!/usr/bin/env zsh

# make folders and link goldstandard
mkdir {1..5}_{1..5}
for mat in {1..5}; do for sample in {1..5}; do
    ln -sf ../../goldstandard/goldstandard_$mat.mat ${mat}_${sample}/goldstandard.mat
done; done

# simulate
~/cwd/src/simulation/simulate.cmd.jl ?_?/

# infer
~/cwd/src/simulation/inference.cmd.jl ?_?/

# roc
~/cwd/roc.R -t ?_?/WP.mat -m ?_?/WP_infer.mat --tf-true ?_?/WT.mat --tf-marker ?_?/WT_infer.mat


