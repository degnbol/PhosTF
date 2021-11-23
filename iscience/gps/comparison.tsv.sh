cut -f1 <(sed 1d family2orn.tsv) | comm - <(sed 1d kinases.family | tr '-' '_' | sort) > comparison.tsv
