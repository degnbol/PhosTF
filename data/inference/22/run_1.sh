export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../PKTFX.jl infer ../../perturbation/logFC_inner.mat 272 153 --epochs 200 --J ../../perturbation/J_inner.mat --lambda 0.01 --lambdaW 0.0 | tee infer_1.log
