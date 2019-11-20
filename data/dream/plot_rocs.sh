for dir in lambda*; do rocauc.py $dir/?_?/sim/WT.mat $dir/?_?/WT_infer_thres.mat -title $dir -out $dir/WT.pdf -trans 0.5 > $dir/auc_WT.log; done
for dir in lambda*; do rocauc.py $dir/?_?/sim/WP.mat $dir/?_?/WP_infer_thres.mat -title $dir -out $dir/WP.pdf -trans 0.5 > $dir/auc_WP.log; done
