#!/usr/bin/env bash
cat ../../raw/amigo2/activators.tsv | cut -f7,11 | cut -d'|' -f1 | sort | uniq > activators.tsv
r activators=repressors
r repressors=positive_elongation
r positive=negative
cat ../../raw/amigo2/positive_regulators.tsv | cut -f3,7,11 | cut -d'|' -f1 | awk '{print $2,$3,$1}' | sed 's/  / /' | tr ' ' '\t' | cut -f-2 | sort | uniq > positive_regulators.tsv
r positive=negative
cat ../../raw/amigo2/protein_phosphorylation.tsv | cut -f7,11 | cut -d'|' -f1 | grep -v "non-specific" > protein_phosphorylation.tsv; cat ../../raw/amigo2/protein_phosphorylation.tsv | cut -f3,7,11 | cut -d'|' -f1 | grep "non-specific" | awk '{print $2 "\t" $1}' >> protein_phosphorylation.tsv; sort protein_phosphorylation.tsv | uniq > temp && mv temp protein_phosphorylation.tsv
r phosphorylation=kinase
cat ../../raw/amigo2/protein_dephosphorylation.tsv | cut -f7,11 | cut -d'|' -f1 | grep -v "MDP-1" > protein_dephosphorylation.tsv; cat ../../raw/amigo2/protein_dephosphorylation.tsv | cut -f3,7,11 | cut -d'|' -f1 | grep "MDP-1" | awk '{print $2 "\t" $1}' >> protein_dephosphorylation.tsv; sort protein_dephosphorylation.tsv | uniq > temp && mv temp protein_dephosphorylation.tsv
r dephosphorylation=phosphatase
