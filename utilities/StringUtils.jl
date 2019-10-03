#!/usr/bin/env julia
module StringUtils

export indent

"""
Empty lines are not indented.
"""
function indent(text::String)
	lines = split(text, '\n')
	@. lines[lines != ""] = '\t' * lines[lines != ""]
	join(lines, '\n')
end



end;

