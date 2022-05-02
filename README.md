# PhosTF
Simulation and inference of networks of gene regulation with Protein Kinase, Transcription Factor and gene elements. Most functions can be called through simulation.jl and inference.jl.

## REQUIRMENTS
- Julia. Currently on 1.7.1. Problems testing on 1.7.2.

## INSTALL
- Run `./install.sh` which sets command `git root`
- Either run `./install.jl` or `./install_instantiate.jl` which installs Julia 
  dependencies as newest versions or as tested versions, respectively.
- Auto-detect environment at startup by adding `startup.jl` to `~/.julia/config/startup.jl`, e.g. by running `cat startup.jl >> ~/.julia/config/startup.jl`

