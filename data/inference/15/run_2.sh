export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../PKTFX.jl infer ../../perturbation/logFC.mat 272 153 --epochs 500 --lambda 1.0 --WT_prior ../../perturbation/WT_prior.mat --J ../../perturbation/J.mat
