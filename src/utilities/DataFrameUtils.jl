using DataFrames

"""
rearranges, does not convert to Matrix of remove row column.
Assumes first column contains rownames.
"""
function subtable(df::DataFrame, colnames, rownames)
    rowcol = names(df)[1]
    @assert rowcol ∉ colnames "Rownames should be in first col named $rowcol, which should not be listed among colnames."
    rowselect(df, rownames)[:, [rowcol; colnames]]
end

"""
rearranges, does not convert to Matrix of remove row column
"""
function rowselect(df::DataFrame, rownames)
    # make sure there are no replicate entries in rownames with the only function
    rowind = [only(findall(df[!, 1] .== k)) for k in rownames]
    df[rowind, :]
end



