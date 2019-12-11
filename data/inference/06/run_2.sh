export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../PKTFX.jl infer ../../pertubation/logFC.mat 272 153 --epochs 500 --J ../../pertubation/J.mat --WT WT_infer_1.mat --WP WP_infer_1.mat
