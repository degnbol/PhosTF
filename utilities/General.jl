#!/usr/bin/env julia
module General
using Logging

export debug

debug() = global_logger(SimpleLogger(stdout, Debug))

end;
