#!/usr/bin/env zsh
# Name conversion, priority is SGD, uniprot, yeastract, then SGD aliases.
# A small disagreement between SGD website and downloaded file is fixed by:
echo $'YFL057C\tAAD16' | cat - ../SGD/gene2ORF.tsv ../uniprot/gene2ORF.tsv ../yeastract/gene2ORF.tsv ../SGD/alias2ORF.tsv > gene2ORF.tsv

# manually add rare unpublished conversions
echo $'YFL023W\tFYV11\nYER036C\tKRE30\nYDR161W\tTCI1' >> gene2ORF.tsv
# all unpublished conversions found on SGD using google (not SGD search) and are referred to in nomenclature history as having been removed as valid names.
# lastly there was a name CUP12 found nowhere but CUP1-2 is a standard name on SGD, so we assume this conversion holds.
echo $'YHR055C\tCUP12' >> gene2ORF.tsv
