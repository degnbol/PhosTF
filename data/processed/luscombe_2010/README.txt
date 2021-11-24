KO_transpose.tsv is the logFC file from raw copied, renamed and . in header have been fixed to - since there are no ORF names with .
# adding ORF at beginning of file to be consistent with style from other filesets.
echo -n $'ORF\t' | cat - TF_KO.tsv > temp && mv temp TF_KO.tsv
r TF=PK
r KO=KO_pval
r PK=TF
