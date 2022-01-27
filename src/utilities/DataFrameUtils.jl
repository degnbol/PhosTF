using DataFrames

# rearranges, does not convert to Matrix of remove row column
subtable(df::DataFrame, colnames, rownames) = rowselect(df[:, colnames], rownames)

# rearranges, does not convert to Matrix of remove row column
function rowselect(df::DataFrame, rownames)
    # make sure there are no replicate entries in rownames with the only function
    rowind = [only(findall(df[!, 1] .== k)) for k in rownames]
    df[rowind, :]
end



