#!/usr/bin/env julia
module ArgParseUtils

using ArgParse
using NamedTupleTools

export parse_arguments

get_field_groups(settings::ArgParseSettings) = [field.group for field in settings.args_table.fields]
get_field_dest_names(settings::ArgParseSettings) = [field.dest_name for field in settings.args_table.fields]

"Returns (positionals, optionals) as named tuples."
function parse_arguments(settings::ArgParseSettings)
    parsed_args = (;parse_args(settings; as_symbols=true)...)
    positional_keys = get_field_dest_names(settings)[get_field_groups(settings) .== "positional"]
    split(parsed_args, Tuple(Symbol.(positional_keys)))
end

"Run the main function f giving it parsed arguments, which are parsed according to the "
function main(settings::ArgParseSettings, f)
    positionals, optionals = parse_arguments(settings)
    @info("Positionals: $positionals")
    @info("Optionals: $optionals")
    f(positionals...; optionals...)
end


end;
