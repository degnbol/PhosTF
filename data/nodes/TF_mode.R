
# assign mode to TF from Genetic Ontology terms (AmiGO2)

# functions
unwhich = function(which, dim=max(which)) {
    y = array(logical(length(which)), dim=dim)
    y[which] = TRUE
    y
}

# read GO terms
activators = read.table("~/cwd/data/processed/amigo2/activators.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)
repressors = read.table("~/cwd/data/processed/amigo2/repressors.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)
positive_reg = read.table("~/cwd/data/processed/amigo2/positive_regulators.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)
negative_reg = read.table("~/cwd/data/processed/amigo2/negative_regulators.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)
positive_elong = read.table("~/cwd/data/processed/amigo2/positive_elongation.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)
negative_elong = read.table("~/cwd/data/processed/amigo2/negative_elongation.tsv", col.names=c("Evidence", "TF"), stringsAsFactors=F)


evidences = unique(c(activators$Evidence, repressors$Evidence, positive_reg$Evidence, negative_reg$Evidence))
# ordered from wrost to best. We don't use ISS, ISA and IBA since they are computational evidence
# http://wiki.geneontology.org/index.php/Guide_to_GO_Evidence_Codes
evidences = c("IEA", "IC", "HMP", "IMP", "IGI", "IDA", "IPI")

# first insert positive and negative regulators so that activators and repressors will replace positive and negative regulators in case of disagreement

get_modes = function(TFs, activator_table, repressor_table) {
    modes = rep("", length(TFs))
    # start at worst, replace with better evidence types
    for (evidence in evidences) {
        activators_found = unwhich(match(activator_table$TF[activator_table$Evidence == evidence], TFs), dim=length(TFs))
        repressors_found = unwhich(match(repressor_table$TF[activator_table$Evidence == evidence], TFs), dim=length(TFs))
        # if a TF is classified as activator and repressor with the same evidence we cannot use it to make a decision
        modes[activators_found & !repressors_found] = "activator"
        modes[repressors_found & !activators_found] = "repressor"
    }
    modes
}

assign_modes = function(TFs) {

    TFs$Mode = get_modes(TFs$TF, activators, repressors)
    TFs$Mode[TFs$Mode==""] = get_modes(TFs$TF, positive_reg, negative_reg)[TFs$Mode==""]
    TFs$Mode[TFs$Mode==""] = get_modes(TFs$TF, positive_elong, negative_elong)[TFs$Mode==""]
    
    # print counts info
    cat(paste(sum(TFs$Mode == "activator"), sum(TFs$Mode == "repressor"), sum(TFs$Mode == "")), "\n")
    
    TFs
}

main = function(infile, outfile) {
    TFs = read.table(infile, col.names="TF", stringsAsFactors=F)
    TFs = assign_modes(TFs)
    write.table(TFs, outfile, sep="\t", quote=F, row.names=F)
}

setwd("~/cwd/data/nodes")
main("TF.txt", "TF_mode.tsv")
main("TF_putative.txt", "TF_putative_mode.tsv")

