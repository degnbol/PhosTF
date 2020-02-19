#!/usr/bin/env bash
./KP2TF_wilcoxon.R `cat $1 | tr '\n' ' ' | cat - <(echo '')`
