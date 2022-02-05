#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import sys
from os.path import basename

# sample
fnames_true = ["../2-insilicoNetworkSimulation/adjacencies/WP_10_1-rep1.adj"]
fnames_score = ["inferredWeights/WP_10_1-rep1.adj"]
fnames = sys.argv[1:]
fnames_true = fnames[:len(fnames)//2]
fnames_score = fnames[len(fnames)//2:]
assert len(fnames_true) == len(fnames_score), (len(fnames_true), len(fnames_score))

for fname_true, fname_score in zip(fnames_true, fnames_score):
    # assert(basename(fname_true) == basename(fname_score))
    trues = pd.read_table(fname_true, sep="\s", header=0, engine="python")
    scores = pd.read_table(fname_score, sep="\s", header=0, engine="python")

    trues = np.asarray(trues.drop(columns="_"))
    scores = np.asarray(scores.drop(columns="_"))

    trues[trues == '.'] = 0
    trues[trues == '+'] = +1
    trues[trues == '-'] = -1

    trues = np.asarray(abs(trues), dtype=bool)
    scores = abs(scores)

    trues = trues.flatten()
    scores = scores.flatten()

    print(roc_auc_score(trues, scores))

