#!/usr/bin/env julia
using DelimitedFiles

outpath = @__DIR__() * "/adjacencies"
ispath(outpath) || mkpath(outpath)

for size in [10, 100]
    for i in 1:5
        adj = fill('.', size, size)
        
        fname = @__DIR__() * "/raw/DREAM4 in-silico challenge/Size $size/Supplementary information/insilico_size$(size)_$(i)/Goldstandard/insilico_size$(size)_$(i)_goldstandard_signed.tsv"
        open(fname) do io
            for line in eachline(io)
                l = split(strip(line))
                source, target = (parse(Int, lstrip(s, 'G')) for s in l[1:2])
                sign = only(l[3])
                adj[target, source] = sign
            end
        end

        writedlm(outpath * "/goldstandard_$(size)_$(i).adj", adj, ' ')
    end
end

