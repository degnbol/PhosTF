for file in ??/score.txt; do cat $file | paste <(echo $file) - ; done
