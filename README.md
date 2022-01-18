# PhosTF
Simulation and inference of networks of gene regulation with Protein Kinase, Transcription Factor and gene elements. Most functions can be called through simulation.jl and inference.jl.

## REQUIRMENTS
- Julia. Currently on 1.6.4

## CONFIG
- Configure some commands by running config.sh.

## INSTALL
- Julia dependencies can be insalled with `install.jl` either by running
  - `./install.jl main` to install newest versions to main environment or
  - `./install.jl env` to install to the isolated environment at git root. Activate it by setting `JULIA_PROJECT`, e.g. temporarily with `. ./activate.sh` or auto-detect environment at startup, e.g. by adding something like `using Pkg; proj = Base.current_project(); proj === nothing || Pkg.activate(proj)` to `~/.julia/config/startup.jl`

