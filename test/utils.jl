using Test
using ONSAS.Utils
using LinearAlgebra

@testset "Utils functions." begin
    eye(3) == Diagonal([1, 1, 1])
end
