#!/usr/bin/env zsh
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WT_10_*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WP_10_*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WT_100_*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
./rocauc.py {../2-insilicoNetworkSimulation/adjacencies,inferredWeights}/WP_100_*.adj | grep -v '^#' | cat <(echo "AUC") - | mlr --tsv stats1 -a mean -f AUC | sed 1d
