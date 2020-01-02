export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../PKTFX.jl infer ../../perturbation/logFC.mat 272 153 --epochs 200 --J ../../perturbation/J.mat --lambda 0.01 --lambdaW 0.0 --WT_prior ../../perturbation/WT_prior.mat | tee infer_1.log
