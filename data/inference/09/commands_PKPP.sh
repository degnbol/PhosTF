#!/usr/bin/env zsh
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
../../../inference.jl "~/cwd/data/perturbation/logFC_inner.mat" 231 199 --J "~/cwd/data/perturbation/J_inner.mat" --epochs 10 --lambda 0.1 --lambdaW 0 --lambdaWT false --WT-prior "" --WP-prior "~/cwd/data/network/KP_mask.mat" --WT-reg "" --WT "~/cwd/data/network/WT.mat" --WP "WP_infer.mat" --trainWT false --quadquad true --PKPP "~/cwd/data/network/PKPP.txt" | tee infer_PKPP.log