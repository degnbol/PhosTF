#!/usr/bin/env julia
module ReadWrite

using DelimitedFiles
using CSV, DataFrames
using Chain: @chain
import JSON3
import BSON
import JLD
# problems with JLD2, not using it

export load, save
export loaddlm, loadmat, savedlm

# for storing network data.
default_identifier = "data"

EXT2DELIM = Dict(".mat"=>' ', ".txt"=>' ', ".ssv"=>' ', ".adj"=>' ', ".csv"=>',', ".tsv"=>'\t')

function ext2delim(fname)
	_, ext = splitext(fname)
    ext ∈ keys(EXT2DELIM) || error("File format not recognized.")
    EXT2DELIM[ext]
end


"""
Load from file using relevant package depending on file extension.
"""
function load(fname::String)
	ext = splitext(fname)[2]
	if ext == ".bson" load_BSON(fname)
	elseif ext == ".jld" load_JLD(fname)
	else loaddlm(fname) end
end
function load(fname::String, cast)
	if splitext(fname)[2] == ".json"
		return load_JSON(fname, cast)
	else
	    return load(fname)
	end
end

function loaddlm(fname::String, T::Union{Type,Nothing}=nothing; header::Bool=false)
    delim = ext2delim(fname)
    if header
        df = CSV.read(fname, DataFrame; delim=delim)
        hasRownames = df |> names |> first |> lowercase in ["_", "rownames", "row"]
        matFirstCol = hasRownames ? 2 : 1
        # if all elements are numbers then there is no strings to parse.
        if !all(eltype.(eachcol(df[!, matFirstCol:end])) .<: Real)
            df[!, matFirstCol:end] = parse_matrix(Matrix(df[!, matFirstCol:end]))
        end
        if T !== nothing
            df[!, matFirstCol:end] = convert.(T, df[!, matFirstCol:end])
        end
        df
    else
        mat = readdlm(fname, delim)
        eltype(mat) <: Real ? mat : parse_matrix(mat)
        T === nothing ? mat : convert.(T, mat)
    end
end

"Load dlm where we try to parse as int, and if that fails as float (which is does automatically)."
function loadmat(fname::String)
	try return loaddlm(fname, Int)
	catch; return loaddlm(fname) end
end

"""
Load from JSON format using JSON3.
For this to work it is necessary that a visible constructor is found that matches the elements of the JSON.
This can mean that you will have to make a constructor that is identical to the "new" function in a struct.
- fname: filename of file to read from
- cast: type to cast the read contents to e.g. a custom struct you have made.
To use this you must define what struct type you are reading to JSON3.
This can be do like so:
JSON3.StructType(::Type{ExampleStruct}) = JSON3.struct()
"""
load_JSON(fname::String, cast) = open(fname) do io return JSON3.read(io, cast) end
load_BSON(fname::String) = BSON.load(fname)[default_identifier]
load_JLD(fname::String) = JLD.load(fname, default_identifier)


"""
Save to file using relevant package depending on file extension.
"""
function save(fname::String, x)
	ext = splitext(fname)[2]
	if ext == ".bson" save_BSON(fname, x)
	elseif ext == ".json" save_JSON(fname, x)
	elseif ext == ".jld" save_JLD(fname, x)
	else savedlm(fname, x) end
end

savedlm(o::Base.TTY, x::AbstractMatrix) = writedlm(o, x)
savedlm(fname::String, x::AbstractMatrix; colnames=nothing, rownames=nothing) = begin
    if colnames === nothing
        writedlm(fname, x, ext2delim(fname))
    else
        df = DataFrame(x, colnames)
        rownames === nothing || insertcols!(df, 1, "_"=>rownames)
        CSV.write(fname, df; delim=ext2delim(fname))
    end
end
# a vector will be written as a column vector.
savedlm(fname, x::AbstractVector) = savedlm(fname, reshape(x, :, 1))

save_JSON(fname::String, x) = open(fname, "w") do io JSON3.write(io, x) end
save_BSON(fname::String, x) = BSON.bson(fname, Dict(default_identifier => x))
save_JLD(fname::String, x) = JLD.save(fname, default_identifier, x)

"""
Deal with string matrices containing '.', '+', '-' to represent 0, 1, -1.
The DelimitedFiles.readdlm should read string characters as Matrix{Any}
"""
parse_matrix(mat::Union{Matrix{Any},Matrix{<:AbstractString}}) = begin
    try
        return parse_matrix(only.(mat))
    catch e
        isa(e, ArgumentError) || rethrow(e)
    end
    mat
end
parse_matrix(mat::Matrix{Char}) = begin
    # all elements has to be '.', '+' or '-'
    all((mat .== '+') .| (mat .== '-') .| (mat .== '.')) || return mat
    out = (mat .== '+') .- (mat .== '-')
end
parse_matrix(mat::Union{BitMatrix,Matrix{<:Real}}) = mat
function pretty_matrix(mat::Matrix{<:Real})
    out = fill('.', size(mat))
    out[mat .== +1] .= '+'
    out[mat .== -1] .= '-'
    out
end

end;
