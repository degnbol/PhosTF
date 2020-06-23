#!/usr/bin/env zsh
#mkdir {1..5}_{1..5}
for mat in {1..5}; do for sample in {1..5}; do
    #ln -sf ../../goldstandard/goldstandard_$mat.mat ${mat}_${sample}/goldstandard.mat
    cd ${mat}_${sample}
    # make WT and WP with {-1, 0, 1} for deactvation, no edge and activation
    #~/cwd/weight.jl random goldstandard.mat 30
    # make net.bson with a fully defined network instance with all it's constants, etc.
    ~/cwd/simulation.jl network
    cd -
done; done
