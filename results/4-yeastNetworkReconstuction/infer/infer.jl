#!/usr/bin/env julia
@src "inference/infer"

iORo, logFC, λB, λW = ARGS

suff = "$iORo$logFC"
infer("logFC_$suff.csv", "TF.txt", "KP.txt", "WT_infer-$suff-$λB-$λW.tsv", "WP_infer-$suff-$λB-$λW.tsv"; 
      delim_mut='_', col_match=r"^\w+", lambda_Bstar=parse(Int, λB), lambda_absW=parse(Int, λW), epochs=10)

