#!/usr/bin/env zsh
sed 1d gky1063_supplemental_files/gps_idmap.txt | mlr --tsv grep 'Saccharomyces cerevisiae' then rename 'Ensembl Gene ID,OLN' then uniq -f 'Family,OLN' | table_unjag.sh 2 $'\t' ';' |
    mlr --tsv put '$Family = sub($Family, ".*/", "")' | table_3D.R Family | mlr --tsv sort -f Family > family2oln.tsv
