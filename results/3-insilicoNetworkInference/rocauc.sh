#!/usr/bin/env zsh
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WT*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WP*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
