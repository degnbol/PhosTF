#!/usr/bin/env julia
using Base.Threads
println("#threads = ", nthreads())

iORo, logFC, λB, λW, WT, WT_mask = [v for ARG in ARGS for v in split(ARG, '+')]
args = join(ARGS, ' ')


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

@src "inference/infer"
infer("../logFC_$suff.csv",
"../TF.txt", "../KP.txt",
"W_infer/WT_infer-$args.tsv", "W_infer/WP_infer-$args.tsv"; 
delim_mut='_', col_match=r"^[\w-]+",
lambda_Bstar=parse(Float64, λB),
lambda_absW=parse(Float64, λW),
epochs=200,
WT=WT, train_WT=WT===nothing,
WT_mask=WT_mask)

