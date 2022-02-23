#!/usr/bin/env zsh
cat aucs{1,2}.tsv > aucs.tsv
mlr --tsv --from aucs.tsv filter '$AUC_T == 1.0 && $AUC_P == 1.0' then cut -x -f AUC_T,AUC_P then uniq -a > good.tsv
mlr --tsv --from aucs.tsv filter '$AUC_T != 1.0 || $AUC_P != 1.0' then cut -x -f AUC_T,AUC_P then uniq -a > bad.tsv
mlr --tsv --from bad.tsv join -j cancel,mean_k,vec,norm,lam,λB,λW --np --ul -f good.tsv > onlygood.tsv
# how resiliant the sim is to different values of reg
mlr --tsv --from onlygood.tsv stats1 -a sum -f AUC_T,AUC_P -g cancel,mean_k,vec,norm,lam then sort -nr AUC_P_sum,AUC_T_sum > reg_resiliance.tsv
