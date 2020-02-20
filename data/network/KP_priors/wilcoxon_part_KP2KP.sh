#!/usr/bin/env bash
./KP2KP_wilcoxon.R `cat $1 | tr '\n' ' ' | cat - <(echo '')`
