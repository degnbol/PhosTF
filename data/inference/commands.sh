# correct annoying line feed
cat hyperparams.tsv | tr -d '\r' > temp && mv temp hyperparams.tsv
sed 1d hyperparams.tsv | tr '\t' ',' | while IFS=',' read folder X J epochs lambda lambdaW lambdaWT WT_mask WP_mask WT_reg WT WP trainWT linex; do
mkdir $folder
echo "#!/usr/bin/env bash
export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../PKTFX.jl infer \"$X\" 272 153 --J \"$J\" --epochs $epochs --lambda $lambda --lambdaW $lambdaW --lambdaWT $lambdaWT --WT_prior \"$WT_mask\" --WP_prior \"$WP_mask\" --WT_reg \"$WT_reg\" --WT \"$WT\" --WP \"$WP\" --trainWT $trainWT --linex $linex | tee infer_1.log" > $folder/commands.sh
chmod 755 $folder/commands.sh
done
