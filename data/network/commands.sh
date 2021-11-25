#!/usr/bin/env zsh
./TF_priors/putative_TF_edges.R
# also makes TF_priors/TF_edges.tsv
./nodes.R
../perturbation/perturbation.R
./TF_mode.R
./WT.R
# make a mask for WT edges with FDR 20
sed 's/-1/-/g' WT_FDR20_sign.mat | sed 's/1/+/g' > WT_FDR20_sign_mask.mat
sed 's/-//g' WT_FDR20_sign.mat > WT_FDR20_abssign_mask.mat

./KP_mask.R
sed 1d KP_mask.csv | cut -f-2 -d, | grep '0$' | sed 's/,0//' > KP_targets_noknownsite.txt
gzip --keep KP_targets_noknownsite.txt

# now to go KP_priors
./KP_priors/commands.sh

# you can update gzipped versions with
ls *.gz | sed 's/\.gz//' | while read line; do gzip -f --keep $line; done

