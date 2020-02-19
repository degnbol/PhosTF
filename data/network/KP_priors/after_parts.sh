grep -hv '^TF' KP2TF_parts/*.tsv | cat <(head -n1 KP2TF_parts/YAL009W.tsv) - > KP2TF_medianest.tsv
# cleanup folder a bit, we can always make this files again with the split_KPs.sh script
rm KP_??.txt
