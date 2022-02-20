#!/usr/bin/env julia
@src "inference/infer"

iORo, logFC, λB, λW, WT, WT_mask = [v for ARG in ARGS for v in split(ARG, '+')]


suff = "$iORo$logFC"
if WT == ""
    WT = nothing
else
    WT = "../WT_$(WT).ssv"
end
if WT_mask == ""
    WT_mask = nothing
else
    WT_mask = "../WT_$(WT_mask).ssv"
end
infer("../logFC_$suff.csv", "../TF.txt", "../KP.txt", "WT_infer-$suff-$λB-$λW.tsv", "WP_infer-$suff-$λB-$λW.tsv"; 
delim_mut='_', col_match=r"^[\w-]+",
lambda_Bstar=parse(Float64, λB),
lambda_absW=parse(Float64, λW),
epochs=200,
WT=WT, train_WT=WT===nothing,
WT_mask=WT_mask)

