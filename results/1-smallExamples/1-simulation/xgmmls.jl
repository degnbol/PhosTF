#!/usr/bin/env julia
@src "visualization/W2graph"

xgmml("net.bson"; o="net.xgmml")
xgmml("net.bson"; o="steady.xgmml", X="steady.tsv")
for mut_id in 1:5
    xgmml("net.bson"; o="steady_mut$mut_id.xgmml", X="steady_mut$mut_id.tsv", highlight=mut_id)
end
xgmml("net.bson"; o="ko.xgmml", X="sim_logFC.tsv")
xgmml("WT_infer.tsv", "WP_infer.tsv"; o="infer.xgmml")


