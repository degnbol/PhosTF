#!/usr/bin/env zsh
nTF=$(cat ../network/TF.txt | wc -l | xargs)
nKP=$(cat ../network/KP.txt | wc -l | xargs)
cat hyperparams.tsv | tr -d '\r' > temp && mv temp hyperparams.tsv
sed 1d hyperparams.tsv | tr '\t' ',' | while IFS=',' read `head -n1 hyperparams.tsv`; do
mkdir $folder
echo "#!/usr/bin/env zsh
export JULIA_NUM_THREADS=12
export OPENBLAS_NUM_THREADS=12
../../../inference.jl \"$X\" $nTF $nKP --J \"$J\" --epochs $epochs --lambda $lambda --lambdaW $lambdaW --lambdaWT $lambdaWT --WT-prior \"$WT_mask\" --WP-prior \"$WP_mask\" --WT-reg \"$WT_reg\" --WT \"$WT\" --WP \"$WP\" --trainWT $trainWT --quadquad $quadquad | tee infer_1.log" > $folder/commands.sh
chmod 755 $folder/commands.sh
done
