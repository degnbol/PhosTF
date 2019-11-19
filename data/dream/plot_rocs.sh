for dir in L*; do rocauc.py $dir/?_1/WT.mat $dir/?_1/WT_infer_thres.mat -title $dir -out $dir/WT.pdf > $dir/auc.log; done
for dir in L*; do rocauc.py $dir/?_1/WP.mat $dir/?_1/WP_infer_thres.mat -title $dir -out $dir/WP.pdf > $dir/auc.log; done
