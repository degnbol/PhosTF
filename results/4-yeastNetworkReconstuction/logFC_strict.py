#!/usr/bin/env python3
import numpy as np
import pandas as pd

# drop columns where a KO is not below zero or a OE is not above.

def drop_signerror(df):
    df = df.copy()
    for colname in df:
        col = df[colname]
        # .1 etc. appended to some columns to make them unique
        orfs = colname.split('.')[0]
        if orfs.endswith(" OE"):
            OE = 1
            orfs = orfs[:-len(" OE")]
        else:
            OE = -1
        orfs = orfs.split("_")
        if not all(col[orfs] * OE > 0):
            print(f"Dropping\n{col[orfs]}\n")
            df.drop(columns=colname, inplace=True)
    return df


dfi = pd.read_csv("logFC_inner.csv", index_col=0)
dfo = pd.read_csv("logFC_outer.csv", index_col=0)

dfid = drop_signerror(dfi)
dfid.to_csv("logFC_inner-strict.csv")
dfod = drop_signerror(dfo)
dfod.to_csv("logFC_outer-strict.csv")
print(f"Removed {dfi.shape[1] - dfid.shape[1]} columns from inner")
print(f"Removed {dfo.shape[1] - dfod.shape[1]} columns from outer")


