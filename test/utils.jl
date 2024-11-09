using Test
using ONSAS.Utils
using LinearAlgebra

@testset "Utils functions." begin
    @test eye(3) == Diagonal([1, 1, 1])
    @test voigt([1 2 3;
                 4 5 6;
                 7 8 9]) == [1, 5, 9, 4, 7, 1, 7, 1, 4]
end
