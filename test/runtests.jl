
using Test
using ONSAS
using SafeTestsets

import ONSAS: nodes2dofs


@safetestset "ONSAS.Materials" begin
    include("materials.jl")
end

@safetestset "ONSAS.BoundaryConditions" begin
    include("boundary_conditions.jl")
end

@safetestset "ONSAS.Mesh" begin
    include("mesh.jl")
end



#using LinearAlgebra
@testset "ONSAS.jl" begin

    @testset "Auxiliar functions" begin
        test_nodes = [1, 5, 3]
        test_dofs = nodes2dofs(test_nodes, 3)
        test_sol_dofs = [1, 2, 3, 13, 14, 15, 7, 8, 9]
        @test isequal(test_dofs, test_sol_dofs) # checks if nodes2dofs is ok
    end

    # @testset "Test: linear von mises problem" begin
    #     include( joinpath( "..", "examples", "vonMisesTruss", "vonMisesTruss.jl" ) )

    #     L0 = sqrt( d^2 + h^2 )

    #     cosalpha = d / L0
    #     sinalpha = h / L0

    #     Ldef12 = norm( [ d+UG[7],h+UG[11]] - [0,0]  )
    #     Ldef23 = norm( [2d,0] - [ d+UG[7],h+UG[11]] )

    #     F12Analy =  +Fx/(2cosalpha) + Fy/(2sinalpha)
    #     F23Analy =  -Fx/(2cosalpha) + Fy/(2sinalpha)

    #     delta12Num = Ldef12-L0 ;
    #     delta23Num = Ldef23-L0 ;

    #     delta12Analy = F12Analy*L0 / (E*A) ;
    #     delta23Analy = F23Analy*L0 / (E*A) ;

    #     @test maximum( abs.( [ delta12Num-delta12Analy, delta12Num-delta12Analy ] ) ) < 1e-8
    # end

end
