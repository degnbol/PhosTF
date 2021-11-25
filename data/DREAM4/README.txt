https://dx.doi.org/10.1371%2Fjournal.pone.0013397

"Gold standard" signed adjacency matrices from the DREAM4 challenge.
In the DREAM4 challenge participants were given in-silico expression levels generated with GeneNetWeaver that was given a signed set of transcription factor interactions.
The goal was to reverse engineer the gold standard set of directed interactions, ignoring sign.
The signed and directed interactions are found in the suppl. under raw/ and reshaped into adjacency matrices in adjacencies/.

An adjacency matrix M has the source node along columns, and the target node along rows, such that (i,j) is from source j to target i. 
This fits naturally with the notion of matrix multiplication, e.g. M^2 shows the targets of source nodes 1 extra step away.
