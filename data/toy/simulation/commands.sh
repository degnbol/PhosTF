./simulation.cmd.jl
for mutnum in {1..4}; do
~/cwd/src/plot_simulation.R 3 2 1 $mutnum
done
~/cwd/inference.jl X_sim.mat 2 3 --epochs 30000
./cytoscape.cmd.jl
