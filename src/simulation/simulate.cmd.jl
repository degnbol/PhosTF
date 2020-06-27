#!/usr/bin/env julia

"""
This command script can be used directly for running multiple simulations or as inspiration.
Folder structure: 
put folders for each run in a working directory. 
Have a goldstandard matrix in each folder. Space delimited table. No header. 0 = no edge, 1 = edge.
Run script from workdir with the folders given as ARGS, 
e.g. if all subfolders are to be run: 
~/cwd/src/simulation/simulate.cmd.jl */
"""

# settings:
max_nₚ = 30
max_attempts = 3  # if a random net turns out to have unstable differential equations.

using Test  # @test_logs
# suppress output of Fire function definitions with ";"
expanduser("~/cwd/weight.jl") |> include;  # random(...) and correct()
expanduser("~/cwd/simulation.jl") |> include;  # other functions

tryrm(fname) = try rm(fname) catch IOError end

for dir in ARGS
    println(dir); cd(dir)
    for attempt in 1:max_attempts
        attempt == 1 || println("attempt $attempt")
        # make WT and WP with {-1, 0, 1} for deactvation, no edge and activation
        random("goldstandard.mat", max_nₚ)
        # there should be nothing to correct
        correct()
        # make net.bson with a fully defined network instance with all it's constants, etc.
        network()

        tryrm("X_sim.mat")
        # simulate logFC values and test that there are no warnings
        try @test_logs (:info, "logFC values simulated") logFC(o="X_sim.mat")
        catch; continue end
        break
    end
    cd("..")
end

