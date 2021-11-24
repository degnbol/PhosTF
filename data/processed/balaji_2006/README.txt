# insert - to correct systematic names
sed 1d ../../raw/balaji_2006/tnet.txt | tr -d '\r' | sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])/\1-\2/' | awk '{print $0 "\t<0.001"}' | cat <(echo $'TF\tTarget\tPval') - > TF_edges.tsv
