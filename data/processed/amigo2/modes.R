
# assign mode of regulation to TFs and KPs from Genetic Ontology terms (AmiGO2)

# functions
unwhich = function(which, dim=max(which)) {
    y = array(logical(dim), dim=dim)
    y[which] = TRUE
    y
}

# read GO terms
setwd("~/cwd/data/processed/amigo2")
activators = read.table("activators.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
repressors = read.table("repressors.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
positive_reg = read.table("positive_regulators.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
negative_reg = read.table("negative_regulators.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
positive_elong = read.table("positive_elongation.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
negative_elong = read.table("negative_elongation.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)

protein_kinase = read.table("protein_kinase.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
protein_phosphatase = read.table("protein_phosphatase.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
protein_phosphorylation = read.table("protein_phosphorylation.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)
protein_dephosphorylation = read.table("protein_dephosphorylation.tsv", col.names=c("Evidence", "Protein"), stringsAsFactors=F)

all_evidences = unique(c(activators$Evidence, repressors$Evidence,
                         positive_reg$Evidence, negative_reg$Evidence,
                         positive_elong$Evidence, negative_elong$Evidence,
                         protein_kinase$Evidence, protein_phosphatase$Evidence,
                         protein_phosphorylation$Evidence, protein_dephosphorylation$Evidence))
# ordered from worst to best
# http://wiki.geneontology.org/index.php/Guide_to_GO_Evidence_Codes
evidences = c("IC", "HMP", "IMP", "TAS", "HDA", "IGI", "IDA", "IPI")
computational_evidence = c("IEA", "IBA", "ISA", "ISS", "ISM") # We don't use IBA, ISA, ISS, and ISM since they are computational evidence
stopifnot(all(all_evidences[!(all_evidences %in% evidences)] %in% computational_evidence))
# remove entries from computational evidence
filter_computational = function(x) {x[!(x$Evidence %in% computational_evidence),]}
activators = filter_computational(activators)
repressors = filter_computational(repressors)
positive_reg = filter_computational(positive_reg)
negative_reg = filter_computational(negative_reg)
positive_elong = filter_computational(positive_elong)
negative_elong = filter_computational(negative_elong)
protein_kinase = filter_computational(protein_kinase)
protein_phosphatase = filter_computational(protein_phosphatase)
protein_phosphorylation = filter_computational(protein_phosphorylation)
protein_dephosphorylation = filter_computational(protein_dephosphorylation)

# first insert positive and negative regulators so that activators and repressors will replace positive and negative regulators in case of disagreement

get_modes = function(proteins, activator_table, repressor_table, activator_string, repressor_string) {
    if (missing(activator_string)) {activator_string="activator"}
    if (missing(repressor_string)) {repressor_string="repressor"}
    modes = rep("", length(proteins))
    # start at worst, replace with better evidence types
    for (evidence in evidences) {
        activators_found = unwhich(match(activator_table$Protein[activator_table$Evidence == evidence], proteins), dim=length(proteins))
        repressors_found = unwhich(match(repressor_table$Protein[repressor_table$Evidence == evidence], proteins), dim=length(proteins))
        # if a protein is classified as activator and repressor with the same evidence we cannot use it to make a decision
        modes[activators_found & !repressors_found | is.na(repressors_found)] = activator_string
        modes[repressors_found & !activators_found | is.na(activators_found)] = repressor_string
    }
    modes
}

get_TF_modes = function(proteins) {
    
    modes = get_modes(proteins, activators, repressors)
    modes[modes==""] = get_modes(proteins, positive_reg, negative_reg)[modes==""]
    modes[modes==""] = get_modes(proteins, positive_elong, negative_elong)[modes==""]
    
    # print counts info
    cat(paste(sum(modes == "activator"), sum(modes == "repressor"), sum(modes == "")), "\n")
    
    modes
}

get_KP_modes = function(proteins) {
    
    modes = get_modes(proteins, protein_kinase, protein_phosphatase, "kinase", "phosphatase")
    modes[modes==""] = get_modes(proteins, protein_phosphorylation, protein_dephosphorylation, "kinase", "phosphatase")[modes==""]
    
    # print counts info
    cat(paste(sum(modes == "kinase"), sum(modes == "phosphatase"), sum(modes == "")), "\n")
    
    modes
}


TFs = unique(c(activators$Protein, repressors$Protein, positive_reg$Protein, negative_reg$Protein, positive_elong$Protein, negative_elong$Protein))
KPs = unique(c(protein_phosphorylation$Protein, protein_dephosphorylation$Protein))
TFs = data.frame(TF=TFs, Mode=get_TF_modes(TFs))
KPs = data.frame(KP=KPs, Mode=get_KP_modes(KPs))


write.table(TFs, "TF2.tsv", sep="\t", quote=F, row.names=F)
write.table(KPs, "KP2.tsv", sep="\t", quote=F, row.names=F)




