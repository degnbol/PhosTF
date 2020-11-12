grep OrderedLocusName ../../raw/uniprot/YEAST_559292_idmapping.dat | cut -f1,3 | sort -u | cat <(echo $'ACC\tOLN') - > acc_oln.tsv
