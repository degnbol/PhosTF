Data downloaded from http://www.yeastract.com/formregmatrix.php where regulatory matrices can be generated.
Each file downloaded with a unique setting. All ORF/Genes selected and all TF from S. Cerevisiae used.
File downloaded is the so-called tsv file, although misleading.
"binding" in filenames refer to binding evidence selected, "expression" refers to expression evidence selected.
"activator" and "inhibitor" refers to selection for TF with respect to expression data.
"and" refers to intersection set. "or" refers to union set. 
"plus" refers to the setting on the website with the same name, that I assumed would be a union as well, but there is 20% extra interactions in the file with "or".
Other overlaps etc. have been tested to make sense, e.g. combining activator and inhibitor files produces the equivalent file that was retrieved without filtering of interaction sign.
Combining files only with binding evidence and only with expression evidence produces the files with either, as expected. 
All files have been preprocessed with a simple cat | sort | uniq

gene2ORF.tsv was created with yeastract name conversion service, although conversion from MALS, MALR, STA1 were manually removed to match the naming of SGD and that in the other files. The only other disagreement with the conversion from ../name_conversion was what AAD16 should be converted to. SGD website says YFL057C but the downloaded from does not. yeastract says YFL057C as well so it is used.

# for goldstandard
comm -12 binding_and_expression.ssv activator_binding_or_expression.ssv > activator_binding_and_expression.ssv
comm -12 binding_and_expression.ssv inhibitor_binding_or_expression.ssv > inhibitor_binding_and_expression.ssv
comm -12 {activator,inhibitor}_binding_and_expression.ssv > unsigned_binding_and_expression.ssv
comm -23 {activator,inhibitor}_binding_and_expression.ssv > only_activator_binding_and_expression.ssv
comm -13 {activator,inhibitor}_binding_and_expression.ssv > only_inhibitor_binding_and_expression.ssv
cat <(sed 's/$/;activator/' only_activator_binding_and_expression.ssv) <(sed 's/$/;inhibitor/' only_inhibitor_binding_and_expression.ssv) <(sed 's/$/;-/' unsigned_binding_and_expression.ssv) | tr ';' '\t' > TF_edges_standard.tsv
../../../src/gene2ORF.R TF_edges_standard.tsv TF_edges.tsv gene2ORF.tsv
echo "TF\tTarget\tRelationship" | cat - TF_edges.tsv > temp && mv temp TF_edges.tsv

# for evaluation of yeast inference
cat binding.ssv | tr ';' '\t' > binding.tsv; ../../../src/gene2ORF.R binding.tsv binding_ORF.tsv ../../name_conversion/gene2ORF.tsv; cat binding_ORF.tsv | sort | uniq > temp && mv temp binding_ORF.tsv
r binding=activator_expression
r activator=inhibitor
comm -12 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > unsigned_expression_ORF.tsv
comm -23 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > only_activator_expression_ORF.tsv
comm -13 <(cat activator_expression_ORF.tsv) <(cat inhibitor_expression_ORF.tsv) > only_inhibitor_expression_ORF.tsv
cat <(sed $'s/$/\tactivator/' only_activator_expression_ORF.tsv) <(sed $'s/$/\tinhibitor/' only_inhibitor_expression_ORF.tsv) <(sed $'s/$/\tambiguous/' unsigned_expression_ORF.tsv) > expression_ORF.tsv
