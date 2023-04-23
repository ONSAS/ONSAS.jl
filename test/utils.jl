using Test: @testset, @test
using ONSAS.Utils

# Test case for ScalarWrapper
@testset "ONSAS.Utils.ScalarWrapper" begin
    # Test initialization of ScalarWrapper
    s = ScalarWrapper(10)
    @test s.x == 10

    # Test getindex method
    @test s[] == 10

    # Test setindex! method
    s[] = 20
    @test s[] == 20
end
