#!/usr/bin/env zsh
./rocauc.py ../2-insilicoNetworkSimulation/adjacencies/WT*.adj inferredWeights/WT*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC
