#!/usr/bin/env python
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import sys
from os.path import isfile
from itertools import product

with open("aucs.tsv", 'w') as fh:
    fh.write("n\ti\trep\tcancel\tmean_k\tvec\tnorm\tBstar\tabsW\tsource\ttarget\tAUC\n")
    
    bool = ["false", "true"]
    regs = ["0.0", "0.1", "0.5", "1.0"]
    for n, i, rep, _cancel, _mean_k, _vec, _norm, Bstar, absW in product([10, 100], range(1, 6), range(1, 6), bool, bool, bool, range(3), regs, regs):
        hyper = f"{_cancel}-{_mean_k}-{_vec}-{_norm}"
        
        for s in ["T", "P"]:
            fname_true = f"../2-insilicoNetworkSimulation/adjacencies/W{s}_{n}_{i}-rep{rep}.adj"
            fname_score = f"inferredWeights/{hyper}/W{s}_{n}_{i}-rep{rep}-{Bstar}-{absW}.adj"
            if not isfile(fname_score):
                print(f'"{fname_score}" not found')
                continue

            trues = pd.read_table(fname_true, sep="\s", header=0, engine="python")
            scores = pd.read_table(fname_score, sep="\s", header=0, engine="python")

            gene_names = list(scores["_"])
            indT = np.char.startswith(gene_names, "T")
            indP = np.char.startswith(gene_names, "P")
            indO = np.char.startswith(gene_names, "O")
            
            trues = np.asarray(trues.drop(columns="_"))
            scores = np.asarray(scores.drop(columns="_"))

            trues[trues == '.'] = 0
            trues[trues == '+'] = +1
            trues[trues == '-'] = -1

            # a dtype=bool conversion fails for some reason
            trues = abs(trues) == 1
            scores = abs(scores)

            for ind, target in zip([indT, indP, indO, indT | indP | indO], ["T", "P", "O", "V"]):
                if not np.any(ind): continue

                try:
                    auc = roc_auc_score(trues[ind, :].flatten(), scores[ind, :].flatten())
                except:
                    # if there is only one class present TODO: there shouldn't be??
                    # indexerror from having dim 42 != 41, idk
                    pass 
                else:
                    fh.write("\t".join([str(v) for v in [n, i, rep, _cancel, _mean_k, _vec, _norm, Bstar, absW, s, target, auc]]) + "\n")


            
