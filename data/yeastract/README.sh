#!/usr/bin/env zsh
# for goldstandard
comm -12 binding_and_expression.ssv activator_binding_or_expression.ssv > activator_binding_and_expression.ssv
comm -12 binding_and_expression.ssv inhibitor_binding_or_expression.ssv > inhibitor_binding_and_expression.ssv
comm -12 {activator,inhibitor}_binding_and_expression.ssv > unsigned_binding_and_expression.ssv
comm -23 {activator,inhibitor}_binding_and_expression.ssv > only_activator_binding_and_expression.ssv
comm -13 {activator,inhibitor}_binding_and_expression.ssv > only_inhibitor_binding_and_expression.ssv
cat <(sed 's/$/;activator/' only_activator_binding_and_expression.ssv) <(sed 's/$/;inhibitor/' only_inhibitor_binding_and_expression.ssv) <(sed 's/$/;-/' unsigned_binding_and_expression.ssv) | tr ';' '\t' > TF_edges_standard.tsv
../gene2ORF/gene2ORF.R TF_edges_standard.tsv TF_edges.tsv gene2ORF.tsv
echo "TF\tTarget\tRelationship" | cat - TF_edges.tsv > temp && mv temp TF_edges.tsv

# for evaluation of yeast inference
cat binding.ssv | tr ';' '\t' > binding.tsv; ../gene2ORF/gene2ORF.R binding.tsv binding_ORF.tsv ../gene2ORF/gene2ORF.tsv; cat binding_ORF.tsv | sort | uniq > temp && mv temp binding_ORF.tsv
r binding=activator_expression
r activator=inhibitor
comm -12 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > unsigned_expression_ORF.tsv
comm -23 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > only_activator_expression_ORF.tsv
comm -13 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > only_inhibitor_expression_ORF.tsv
cat <(sed $'s/$/\tactivator/' only_activator_expression_ORF.tsv) <(sed $'s/$/\tinhibitor/' only_inhibitor_expression_ORF.tsv) <(sed $'s/$/\tambiguous/' unsigned_expression_ORF.tsv) > expression_ORF.tsv
