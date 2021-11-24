The dataset comprised of processed data produced in Zelezniak at al, Cell Systems 2018 study. 

metabolites_dataset.data_prep.tsv - processed, batch corrected metabolite levels from glycolysis, pentose phosphate and TCA cycle pathways of yeast kinase knockouts
columns:
  metabolite_id - metabolite ID
  kegg_id - corresponding KEGG metabolite ID (Kanehisa et al 2016)  
  official_name - the common name of metabolite
  dataset - the protocol that was used to generate the dataset, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018 
  KO_ORF - kinase knockout's ORF, where WT - parental strain/wild-type, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018
  KO_gene - kinase knockout's yeast gene name, where WT - parental strain/wild-type, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018
  replicate - biological replicate
  value - metabolite signal obtained from SRM-MS/MS experiment, corrected for batch effects, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018

proteins_dataset.data_prep.tsv - processed, batch corrected protein levels of yeast kinase knockouts
columns:
  ORF - ORF of measured protein
  sample - sample name from the study
  replicate - replicate
  KO_ORF - kinase knockout's ORF, where WT - parental strain/wild-type, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018
  KO_gene - kinase knockout's yeast gene name, where WT - parental strain/wild-type, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018
  value - label-free protein signal quantification using SWATH-MS, corrected for batch effects, for more details please see STAR methods in Zelezniak et al, Cell Systems, 2018

References:
  Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.; KEGG as a reference resource for gene and protein annotation. Nucleic Acids Res. 44, D457-D462 (2016)
