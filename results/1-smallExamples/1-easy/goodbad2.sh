#!/usr/bin/env zsh
mlr --tsv --from aucs3.tsv put '$AUC_sum = $AUC_P + $AUC_T + $AUC_Ppos + $AUC_Tpos + $AUC_Pneg + $AUC_Tneg' > aucs_sums.tsv
mlr --tsv --from aucs_sums.tsv filter '$AUC_sum == 6' then cut -x -r -f 'AUC_*' then uniq -a > good.tsv
mlr --tsv --from aucs_sums.tsv filter '$AUC_sum != 6' then cut -x -r -f 'AUC_*' then uniq -a > bad.tsv
mlr --tsv --from bad.tsv join -j cancel,mean_k,vec,norm,lam,λB,λW --np --ul -f good.tsv > onlygood.tsv
# how resiliant the sim is to different values of reg
mlr --tsv --from aucs_sums.tsv stats1 -a sum --fr 'AUC_*' -g cancel,mean_k,vec,norm,lam then sort -nr AUC_sum_sum,AUC_P_sum,AUC_T_sum,AUC_Ppos_sum,AUC_Tpos_sum,AUC_Pneg_sum,AUC_Tneg_sum > reg_resiliance.tsv
mlr --tsv --from aucs_sums.tsv stats1 -a sum --fr 'AUC_*' -g cancel,mean_k,vec,norm,lam then grep -v "error" then \
    sort -nr AUC_sum_sum,AUC_P_sum,AUC_T_sum,AUC_Ppos_sum,AUC_Tpos_sum,AUC_Pneg_sum,AUC_Tneg_sum > AUCsums.tsv

mlr --tsv --from aucs_sums.tsv filter '$AUC_sum == 6' then cut -x -r -f 'AUC_*,λB,λW' then uniq -a -c then sort -nr count > nPerf.tsv
