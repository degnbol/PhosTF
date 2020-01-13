Collect all measured nodes and categorize them as KP, activator TF, repressor TF, or alternatively as non-regulator gene.
Purpose is to allow any node that is regulator to have regulatory effects, even if it is never perturbated.
nodes.R is used first to get the lists of KPs, TFs, and V, where TFs can be given mode from GO terms with TF_modes.R.
perturbation.R in another folder has to be run to make perturbation tables before we can use TF_modes_fill.R to fill in the remaining regulatory modes, since perturbation data is used as fallback to infer mode.
KP.txt is known KPs, KP_KO.txt is KPs mutated in at least 1 perturbation experiment.
Collect all TF edges with binding evidence with provided scores, mainly p-values.
Purpose is to create a prior for KP edge inference, TF edges should not be trained on at all.
First make node sets in nodes folder. Then make unsigned edges with unsigned_edges.R. Then sign the remaining TFs back in nodes folder. Then sign all the edges with signed_edge.R.
