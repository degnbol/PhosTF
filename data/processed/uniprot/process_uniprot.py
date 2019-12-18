#!/usr/bin/env python3
import re

infname = "uniprot_yeast.txt"
outfname = "gene2ORF.tsv"

with open(infname, 'r') as infile:
    conversions = []
    alt_conversions = []

    for line in infile:
        # split by 2 spaces, we want 1 space to be allowed within a cell
        row = line.strip().split("  ")
        cells = []
        for cell in row:
            cell = cell.strip()
            if cell != '': cells.append(cell)
        
        genes = re.split(';|/', cells[0])
        ORFs = re.split(';|/', cells[1])
        conversions.append([ORFs[0], genes[0]])
        for alt_gene in genes[1:]:
            for alt_ORF in ORFs[1:]:
                alt_conversions.append([alt_ORF, alt_gene])

with open(outfname, 'w') as outfile:
    # put main conversions first so the first match found is the best conversion
    for conversion in conversions + alt_conversions:
        outfile.write(conversion[0].strip() + "\t" + conversion[1].strip() + "\n")
