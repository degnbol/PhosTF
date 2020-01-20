#!/usr/bin/env julia
module CLI

export inout
export stdread, stdwrite

"""
Get (i, o) from two positional arguments: io, o.
They should default to nothing.

Example:
```
@main function main(io=nothing, o=nothing)
	i, o = inout(io, o)
	...
end
```
"""
inout(io::Nothing, o::Nothing) = stdin, stdout
inout(io::String, o::String) = io, o
inout(io::String, o::Nothing) = _inout(stdin, io)
_inout(i::Base.PipeEndpoint, o::String) = i, o
_inout(_::Base.TTY, i::String) = i, stdout



function stdread(i)
	if i == stdin return read(i, String)
	else open(i) do io
		return read(io, String)
	end end
end

function stdwrite(o, x)
	if o == stdout return write(o, x)
	else open(o, write=true) do io
		return write(io, x)
	end end
end


end;