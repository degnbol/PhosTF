./TF_priors/putative_TF_edges.R
./nodes.R
../perturbation/perturbation.R
./TF_mode.R
./TF_priors/TF_edge_weights.R
./edge_counts.R
./WT.R
# make PKPP
sed 1d KP.tsv | cut -f2 | sed 's/phosphatase/-1/' | sed 's/kinase/1/' | sed 's/^$/0/' > PKPP.txt
