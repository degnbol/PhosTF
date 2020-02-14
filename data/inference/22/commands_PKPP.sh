#!/usr/bin/env zsh
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
../../../inference.jl "~/cwd/data/perturbation/logFC_outer.mat" 231 199 --J "~/cwd/data/perturbation/J_outer.mat" --epochs 10 --lambda 0.1 --lambdaW 0 --lambdaWT false --WT-prior "~/cwd/data/network/WT_mask.mat" --WP-prior "" --WT-reg "" --WT "WT_infer.mat" --WP "WP_infer.mat" --trainWT true --quadquad false --PKPP "~/cwd/data/network/PKPP.txt" | tee infer_PKPP.log
