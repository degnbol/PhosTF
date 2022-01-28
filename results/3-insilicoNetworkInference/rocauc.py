#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import sys

fnames_true = ["../2-insilicoNetworkSimulation/adjacencies/WP_10_1-rep1.adj"]
fnames_score = ["inferredWeights/WP_10_1-rep1.adj"]
fnames = sys.argv[1:]
fnames_true = fnames[:len(fnames)//2]
fnames_score = fnames[len(fnames)//2:]

for fname_true, fname_score in zip(fnames_true, fnames_score):
    print("# ", fname_true, "\t", fname_score)
    trues = pd.read_table(fname_true, sep=' ', header=0)
    scores = pd.read_table(fname_score, sep=' ', header=0)

    trues = np.asarray(trues.drop(columns="_"))
    scores = np.asarray(scores.drop(columns="_"))

    trues[trues == '.'] = 0
    trues[trues == '+'] = +1
    trues[trues == '-'] = -1
    trues = np.asarray(abs(trues), dtype=bool)

    trues = trues.flatten()
    scores = scores.flatten()

    print(roc_auc_score(trues, scores))

