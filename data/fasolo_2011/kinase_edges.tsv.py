#!/usr/bin/env python3
# REQUIRES: raw/FasoloSuppTableS1.pdf.tsv
# WRITES: kinase_edges.tsv

kinases = []
targets = []

with open("raw/FasoloSuppTableS1.pdf.docx.tsv") as fh:
    # first line contains kinases
    ks = next(fh).strip().split('\t')
    print("kinases:", ks)
    
    for line in fh:
        if line.strip(" \t\r\n") == "":
            # new table, new kinases
            ks = next(fh).strip().split('\t')
            print("kinases:", ks)
        else:
            ts = line.strip().split('\t')
            assert len(ks) >= len(ts), f"len({ks}) < len({ts})"
            kinases.extend(ks[:len(ts)]) # the row may not have as many elements as the header
            targets.extend(ts) # some entries will be ""

with open("kinase_edges.tsv", 'w') as fh:
    fh.write("kinase\ttarget\n")
    for kinase, target in zip(kinases, targets):
        if target != "":
            fh.write(f"{kinase}\t{target}\n")


