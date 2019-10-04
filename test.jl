#!/usr/bin/env julia
module Test

function test01(; kwargs...)
	for (k,v) in kwargs
		println(k, ",", v)
	end
end
	
function test02(; kwargs...)
	test01(;kwargs...)
end

test02(a="davs")


end;


