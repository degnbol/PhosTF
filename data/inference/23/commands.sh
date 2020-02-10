#!/usr/bin/env zsh
export JULIA_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
../../../inference.jl "~/cwd/data/perturbation/logFC_outer.mat" 231 199 --J "~/cwd/data/perturbation/J_outer.mat" --epochs 100 --lambda 0.1 --lambdaW 0 --lambdaWT true --WT-prior "" --WP-prior "" --WT-reg "~/cwd/data/network/WT_p.mat" --WT "~/cwd/data/network/WT.mat" --WP "~/cwd/data/network/KP_priors/WP_noise.mat" --trainWT true --quadquad false | tee infer_1.log
