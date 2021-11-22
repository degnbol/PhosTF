#!/usr/bin/env julia
module UnitTest
using Test
include("ArrayUtils.jl")
using .ArrayUtils

v = rand(3)

@testset "ArrayUtils" begin
	@test hcatpad(v) == hcat(v)
	@test hcatpad([v, v]) == hcat(v, v)
	@test size(reorder(rand(4,4), [3,4,1,2])) == (4,4)
	@test_logs (:warn, "Indexing outside range of shortest axis.") reorder(rand(5,4), [3,4,1,2,5])
	@test_throws BoundsError reorder(rand(4,4), [3,4,1,2,5])
	@test size(reorder(rand(4,4), [3,4,1])) == (3,3)
	@test tological([2, 1], 3) == [1, 1, 0]
	@test tological([2],    3) == [0, 1, 0]
	@test tological( 2 ,    3) == [0, 1, 0]
end

end;