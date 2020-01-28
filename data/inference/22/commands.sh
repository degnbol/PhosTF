#!/usr/bin/env zsh
export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../inference.jl "~/cwd/data/perturbation/logFC_outer.mat" 231 199 --J "~/cwd/data/perturbation/J_outer.mat" --epochs 200 --lambda 0.1 --lambdaW 0 --lambdaWT false --WT-prior "~/cwd/data/network/WT_mask.mat" --WP-prior "" --WT-reg "" --WT "~/cwd/data/network/WT.mat" --WP "~/cwd/data/network/KP_priors/WP.mat" --trainWT true --linex false | tee infer_1.log
