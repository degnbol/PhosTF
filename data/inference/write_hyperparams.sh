#!/usr/bin/env zsh

# Write individual command files for a run on each line in $1 where the file is named $2 and the log is saved in $3. E.g. hyperparams.tsv commands.sh infer.log

cat $1 | tr -d '\r' > temp && mv temp $1

nTF=$(cat ../network/TF.txt | wc -l | xargs)
nKP=$(cat ../network/KP.txt | wc -l | xargs)
sed 1d $1 | tr '\t' ',' | while IFS=',' read `head -n1 $1`; do
mkdir -p $folder
echo "#!/usr/bin/env zsh
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
../../../inference.jl \"$X\" $nTF $nKP --J \"$J\" --epochs $epochs --lambda $lambda --lambdaW $lambdaW --lambdaWT $lambdaWT --WT-prior \"$WT_mask\" --WP-prior \"$WP_mask\" --WT-reg \"$WT_reg\" --WT \"$WT\" --WP \"$WP\" --trainWT $trainWT | tee $3" > $folder/$2
chmod 755 $folder/$2
done
