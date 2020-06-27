#!/usr/bin/env julia
using Test  # @test_logs
# suppress output of Fire function definitions with ";"
expanduser("~/cwd/weight.jl") |> include;  # random(...) and correct()
expanduser("~/cwd/simulation.jl") |> include;  # other functions

tryrm(fname) = try rm(fname) catch IOError end

for mat in 1:5
    for sample in 1:5
        dir="$(mat)_$(sample)"; println(dir); cd(dir)
        for attempt in 1:10
            attempt == 1 || println("attempt $attempt")
            # make WT and WP with {-1, 0, 1} for deactvation, no edge and activation
            random("goldstandard.mat", 30)
            # there should be nothing to correct
            correct()
            # make net.bson with a fully defined network instance with all it's constants, etc.
            network()
            xgmml(; o="net.xgmml")

            # simulate()
            # plot(nₚ, nₜ)  # we could plot now (nₚ and nₜ are not defined, just an example)
            # steadystate()
            
            tryrm("X_sim.mat")
            # simulate logFC values and test that there are no warnings
            try @test_logs (:info, "logFC values simulated") logFC(o="X_sim.mat")
            catch; continue end
            break
        end
        cd("..")
    end
end
