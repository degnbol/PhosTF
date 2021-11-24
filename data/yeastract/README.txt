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

