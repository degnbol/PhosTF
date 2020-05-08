./TF_priors/putative_TF_edges.R
./nodes.R
../perturbation/perturbation.R
./TF_mode.R
./WT.R
# make PKPP
sed 1d KP.tsv | cut -f2 | sed 's/phosphatase/-1/' | sed 's/kinase/1/' | sed 's/^$/0/' > PKPP.txt
# make a mask for WT edges with FDR 20
sed 's/-1/-/g' WT_FDR20_sign.mat | sed 's/1/+/g' > WT_FDR20_sign_mask.mat
sed 's/-//g' WT_FDR20_sign.mat > WT_FDR20_abssign_mask.mat
