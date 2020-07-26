#!/usr/bin/env zsh

# Write individual command files for a run on each line in $1 where the file is named $2 and the log is saved in $3. E.g. hyperparams.tsv commands.sh infer.log

cat $1 | tr -d '\r' > temp && mv temp $1

nTF=$(cat ../network/TF.txt | wc -l | xargs)
nKP=$(sed 1d ../network/KP_protein.tsv | wc -l | xargs)
echo "nTF=$nTF\tnKP=$nKP"

sed 1d $1 | tr '\t' ',' | while IFS=',' read `head -n1 $1`; do
echo "folder: $folder"
mkdir -p $folder
echo "#!/usr/bin/env zsh
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
../../../inference.jl \"$X\" $nTF $nKP --J \"$J\" --epochs $epochs --opt $opt -η $lr -d $decay -λ $lambda --lambdaW $lambdaW --lambdaWT $lambdaWT --WT-prior \"$WT_mask\" --WP-prior \"$WP_mask\" --WT-reg \"$WT_reg\" --WT \"$WT\" --WP \"$WP\" --trainWT $trainWT | tee $3" > $folder/$2
chmod 755 $folder/$2
done
