`adjacencies.jl` was run first to generate adjacency matrices from DREAM4 gold standard matrices.
Here they have 30 percent of TFs repurposed as KPs as opposed to creating new KP->TFKP edges with a similar approach to how the TF gold standard matrices were generated in DREAM4.
`GNWPhosNets.jl` is run next to generate GeneNetWeaverPhos networks with various network constants, decay rates, etc. set randomly.
`sim_logFCs.jl` is run last to simulate logFC values for these insilico networks.
