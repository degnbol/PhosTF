#!/usr/bin/env zsh
mkdir {1..5}_{1..5}
for mat in {1..5}; do for sample in {1..5}; do ln -sf ../../goldstandard/goldstandard_$mat.mat ${mat}_${sample}/goldstandard.mat; done; done
