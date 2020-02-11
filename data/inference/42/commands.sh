#!/usr/bin/env zsh
export JULIA_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
../../../inference.jl "~/cwd/data/perturbation/logFC_inner.mat" 231 199 --J "~/cwd/data/perturbation/J_inner.mat" --epochs 100 --lambda 0.1 --lambdaW 0 --lambdaWT true --WT-prior "~/cwd/data/network/WT_mask.mat" --WP-prior "~/cwd/data/network/KP_mask.mat" --WT-reg "" --WT "~/cwd/data/network/WT.mat" --WP "" --trainWT true --quadquad false | tee infer_1.log
