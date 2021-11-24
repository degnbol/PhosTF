#!/usr/bin/env zsh
# Extract columns with standard and systematic names and remove incomplete entries:
cat raw/SGD_features.tab | cut -d$'\t' -f4,5 | sed -e '/^[[:space:]]$/d' | sed -e 's/[[:space:]]$//' | sed $'/\t/!d' > gene2ORF.tsv
# We can check that all conversions are unique:
wc -l gene2ORF.tsv
cat gene2ORF.tsv | sort | uniq | wc -l
cat gene2ORF.tsv | cut -d$'\t' -f1 | sort | uniq | wc -l
cat gene2ORF.tsv | cut -d$'\t' -f2 | sort | uniq | wc -l
# They all have the same length.

sed 's/ .*//' raw/orf_trans_all.fasta | sed $'/>/s/$/\t/' > orfs.aa
makeblastdb -in ../orfs.aa -dbtype prot
blastp -num_threads 4 -evalue 1e-10 -outfmt 6 -db orfs.aa -query orfs.aa -out orfs_blast.tab
cat orfs_blast.tab | sort -k12nr > orfs_blast_sort.tab

grep -v '^!' raw/gene_association.sgd | cut -f3,11 | cut -d'|' -f1 | grep -v '\t$' | awk '{print $2 "\t" $1}' | sort | uniq > gene_association_gene2ORF.tsv

# collect GO terms (filter out lines where GO term is not indicated)
grep GO raw/go_slim_mapping.tab | cut -f1,4,6 > GO.tsv
