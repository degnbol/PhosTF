#!/usr/bin/env zsh
mkdir -p raw
cd raw
# Reference genome downloaded from https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/
wget 'http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz'
# gene associations downloaded from http://downloads.yeastgenome.org/?prefix=curation/literature/
wget 'http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'
