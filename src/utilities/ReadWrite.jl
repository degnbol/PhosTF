#!/usr/bin/env julia
module ReadWrite

using DelimitedFiles
import JSON3
import BSON
import JLD
using Flux: TrackedArray, Tracker.data
# problems with JLD2, not using it

export load, save
export loaddlm, loadmat, savedlm

default_identifier = "data"

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
	else return load(fname)
	end
end

function loaddlm(fname::String)
	ext = splitext(fname)[2]
	if ext in [".mat", ".txt", ".ssv"] readdlm(fname, ' ')
	elseif ext == ".csv" readdlm(fname, ',')
	elseif ext == ".tsv" readdlm(fname, '\t')
	else error("File format not recognized.") end
end
function loaddlm(fname::String, T::Type)
	ext = splitext(fname)[2]
	if ext in [".mat", ".txt", ".ssv"] readdlm(fname, ' ', T)
	elseif ext == ".csv" readdlm(fname, ',', T)
	elseif ext == ".tsv" readdlm(fname, '\t', T)
	else error("File format not recognized.") end
end
"Load dlm where we try to parse as int, and if that fails as float (which is does automatically)."
function loadmat(fname::String)
	try return loaddlm(fname, Int64)
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

function savedlm(fname::String, x::AbstractArray)
	ext = splitext(fname)[2]
	if ext in [".mat", ".txt", ".ssv"] writedlm(fname, x, ' ')
	elseif ext == ".csv" writedlm(fname, x, ',')
	elseif ext == ".tsv" writedlm(fname, x, '\t')
	else error("File format not recognized.") end
end
savedlm(fname::String, x::TrackedArray) = savedlm(fname, data(x))
savedlm(o::Base.TTY, x::Matrix) = writedlm(o, x)
savedlm(o::Base.TTY, x::TrackedArray) = writedlm(o, data(x))

save_JSON(fname::String, x) = open(fname, "w") do io JSON3.write(io, x) end
save_BSON(fname::String, x) = BSON.bson(fname, Dict(default_identifier => x))
save_JLD(fname::String, x) = JLD.save(fname, default_identifier, x)

end;
