#!/usr/bin/env zsh
# download the necessary files from biogrid
./raw.sh

# only use yeast entries
mlr --tsvlite --from raw/BIOGRID-PTM.ptmtab.txt.gz filter '${Organism Name} == "Saccharomyces cerevisiae (S288c)"' then \
    cut -x -f 'Organism ID,Organism Name' > yeast_PTM.tsv
mlr --tsvlite --from raw/BIOGRID-PTM-RELATIONSHIPS.ptmrel.txt.gz filter '${Organism Name} == "Saccharomyces cerevisiae (S288c)"' then \
    cut -x -f 'Organism ID,Organism Name' > yeast_PTM_relationships.tsv

# make P_edges.tsv and regulatory_edges.tsv
./PTM2edges.R

# extract other parts of the tables
mlr --tsvlite --from yeast_PTM.tsv filter '${Post Translational Modification} == "Phosphorylation"' then \
    uniq -f 'Systematic Name' | sed 1d | sort > targets.txt
mlr --tsvlite --from yeast_PTM.tsv uniq -f "Systematic Name,Position,Residue" then sort -f "Systematic Name" > PTM.res
mlr --tsvlite --from yeast_PTM.tsv uniq -f "Systematic Name,Sequence" then \
    filter '$Sequence != "-" && ${Systematic Name} != "-"' then sort -f "Systematic Name" > sequences.linfa
sed 1d sequences.linfa | sed 's/^/>/' | tr '\t' '\n' > sequences.fa

mlr --tsv --from yeast_PTM_relationships.tsv uniq -f "Systematic Name,Relationship" then sort -f "Systematic Name" > KP.tsv
mlr --tsv --from KPs.tsv filter '$Relationship == "kinase"' | sed 1d | sort > kinases.txt
mlr --tsv --from KPs.tsv filter '$Relationship == "phosphatase"' | sed 1d | sort > phosphatases.txt

echo -n "number of unique sites = "
mlr --tsv --from yeast_PTM_relationships.tsv filter '$Relationship != "-"' then uniq -f '#PTM ID' then count | sed 1d
