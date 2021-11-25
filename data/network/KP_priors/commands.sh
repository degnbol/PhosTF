#!/usr/bin/env zsh
./KP2TF_wilcoxon_setup.R  # makes KP.txt and wilcoxon.RData
./split_KPs.sh

# make or maybe remove it first if already present
mkdir KP2TF_parts
# do in parallel:
./wilcoxon_part.sh KP_01.txt
./wilcoxon_part.sh ...

./after_parts.sh
./KP2TF.R
./WP.R

# make sure to update the zipped versions for work across platforms.
ls *.gz | sed 's/\.gz//' | while read line; do gzip -cf $line > ${line}.gz; done
