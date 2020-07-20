for file in ??/score.txt; do cat $file | paste <(echo $file) - ; done
for file in ??/GO_score.txt; do cat $file | paste <(echo $file) - ; done
