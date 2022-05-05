#!/usr/bin/env julia
# Run me from one of the folders that has the WT.adj and WP.adj
@src "inference/infer"
@src "visualization/W2graph"

println("infer ideal")

for λBstar in [1., .1, .01, .001, 0.]
    for λabsW in [1., .1, .01, .001, 0.]
        WT_infer = "WT_inferIdeal-B$(λBstar)-W$(λabsW).tsv"
        WP_infer = "WP_inferIdeal-B$(λBstar)-W$(λabsW).tsv"
        cy_infer =    "inferIdeal-B$(λBstar)-W$(λabsW).xgmml"
        # infer("../ideal_logFC.tsv", r"T.*", r"P.*", WT_infer, WP_infer; epochs=20000, lambda_Bstar=λBstar, lambda_absW=λabsW)
        xgmml(WT_infer, WP_infer; o=cy_infer, typecol=false, showlab=false)
    end
end


println("infer sims")
@src "simulation/simulate";

WT = "../WT.adj"
WP = "../WP.adj"

for sim in 1:5
    net = "net$(sim).bson"
    network(WT, WP; header=true, o=net)
    logFC(net, "sim$(sim)_logFC.tsv")
    
    # for λBstar in [1., .1, .01, .001, 0.]
    for λBstar in [10., 1., .1, .01, 0.]
        # for λabsW in [1., .1, .01, .001, 0.]
        for λabsW in [0.]
            WT_infer = "WT_infer$(sim)-B$(λBstar)-W$(λabsW).tsv"
            WP_infer = "WP_infer$(sim)-B$(λBstar)-W$(λabsW).tsv"
            cy_infer =    "infer$(sim)-B$(λBstar)-W$(λabsW).xgmml"
            infer("sim$(sim)_logFC.tsv", r"T.*", r"P.*", WT_infer, WP_infer; epochs=20000, lambda_Bstar=λBstar, lambda_absW=λabsW)
            xgmml(WT_infer, WP_infer; o=cy_infer, typecol=false, showlab=false)
        end
    end
end


