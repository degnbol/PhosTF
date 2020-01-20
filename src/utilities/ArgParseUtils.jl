#!/usr/bin/env julia
module ArgParseUtils

using ArgParse
using NamedTupleTools

export parse_arguments

get_field_groups(settings) = [field.group for field in settings.args_table.fields]
get_field_dest_names(settings) = [field.dest_name for field in settings.args_table.fields]

"Returns (positionals, optionals) as named tuples."
function parse_arguments(settings)
    parsed_args = (;parse_args(settings; as_symbols=true)...)
    positional_keys = get_field_dest_names(settings)[get_field_groups(settings) .== "positional"]
    split(parsed_args, Tuple(Symbol.(positional_keys)))
end

function main(settings, f)
    positionals, optionals = parse_arguments(settings)
    @info("Positionals: $positionals")
    @info("Optionals: $optionals")
    f(positionals...; optionals...)
end


end;
