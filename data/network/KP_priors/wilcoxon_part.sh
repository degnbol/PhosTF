#!/usr/bin/env zsh
kpnames=$(cat $1 | tr '\n' ' ' | cat - <(echo ''))
./KP2TF_wilcoxon.R $kpnames
