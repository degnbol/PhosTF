#!/usr/bin/env julia
module General
using Logging

export debug
export abspath_

debug() = global_logger(SimpleLogger(stdout, Debug))

abspath_(path::AbstractString) = abspath(expanduser(path))

end;
