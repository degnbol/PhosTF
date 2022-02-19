#!/usr/bin/env julia
@src "inference/infer"

iORo, logFC, λB, λW = [v for ARG in ARGS for v in split(ARG, '_')]

suff = "$iORo$logFC"
infer("../logFC_$suff.csv", "../TF.txt", "../KP.txt", "WT_infer-$suff-$λB-$λW.tsv", "WP_infer-$suff-$λB-$λW.tsv"; 
delim_mut='_', col_match=r"^[\w-]+", lambda_Bstar=parse(Float64, λB), lambda_absW=parse(Float64, λW), epochs=10)

