From "Diverse protein kinase interactions identified by protein microarrays reveal novel connections between cellular processes" by Fasolo 2011.
Supplementary table S1 downloaded from http://genesdev.cshlp.org/content/suppl/2011/03/30/25.7.767.DC1
Contains 1023 PK interactions.

# correct sys names by inserting a -
sed -E 's/(Y[A-Z]{2}[0-9]{3}[A-Z])([A-Z])/\1-\2/' P_edges_ORF.tsv > temp && mv temp P_edges_ORF.tsv
