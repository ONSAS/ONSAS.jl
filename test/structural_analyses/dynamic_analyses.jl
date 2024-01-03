using Test, Dictionaries, SparseArrays, LinearSolve, LinearAlgebra

using ONSAS.StructuralAnalyses
using ONSAS.DynamicStates
using ONSAS.Assemblers
using ONSAS.StructuralSolvers
using ONSAS.Nodes
using ONSAS.Trusses
using ONSAS.Circles

d = 1.0
h = 2.0
## Nodes
n1 = Node(0.0, 0.0, 0.0)
n2 = Node(d, 0.0, h)
n3 = Node(2d, 0.0, 0.0)
nodes = [n1, n2, n3]
## Entities
truss1 = Truss(n1, n2, Circle(d / 10), "left_truss") # [n1, n2]
truss2 = Truss(n2, n3, Circle(d / 10), "right_truss") # [n2, n3]

# FullDynamicState
free_dofs = [3, 4, 5]
ΔUᵏ = rand(3)
Uᵏ = rand(9)
Udotᵏ = rand(9)
Udotdotᵏ = rand(9)
Fₑₓₜᵏ = rand(9)
Fᵢₙₜᵏ = rand(9)
Fᵢₙₑᵏ = rand(9)
Fᵥᵢₛᵏ = rand(9)
Mᵏ = spzeros(9, 9)
Cᵏ = spzeros(9, 9)
Kᵏ = spzeros(9, 9)
Kₛᵏ = spzeros(9, 9)
res_forces = zeros(3)
ϵᵏ = dictionary([truss1 => rand(3, 3), truss2 => rand(3, 3)])
σᵏ = dictionary([truss1 => rand(3, 3), truss2 => rand(3, 3)])
assembler = Assembler(2)
iter_state = ResidualsIterationStep()
linear_system = init(LinearProblem(Kₛᵏ, res_forces))

sst_rand = FullDynamicState(free_dofs,
                            ΔUᵏ, Uᵏ, Udotᵏ, Udotdotᵏ,
                            Fₑₓₜᵏ, Fᵢₙₜᵏ, Fᵢₙₑᵏ, Fᵥᵢₛᵏ,
                            Kᵏ, Mᵏ, Cᵏ, Kₛᵏ, res_forces,
                            ϵᵏ, σᵏ,
                            assembler,
                            iter_state,
                            linear_system)

sst_rand_sol = DynamicState(Uᵏ, Udotᵏ, Udotdotᵏ, ϵᵏ, σᵏ)

@testset "ONSAS.StructuralAnalyses.DynamicAnalyses.FullDynamicState" begin
    @test displacements(sst_rand) == Uᵏ
    @test velocity(sst_rand) == Udotᵏ
    @test acceleration(sst_rand) == Udotdotᵏ
    @test Δ_displacements(sst_rand) == ΔUᵏ
    @test external_forces(sst_rand) == Fₑₓₜᵏ
    @test internal_forces(sst_rand) == Fᵢₙₜᵏ
    @test viscous_forces(sst_rand) == Fᵥᵢₛᵏ
    @test inertial_forces(sst_rand) == Fᵢₙₑᵏ
    @test tangent_matrix(sst_rand) == Kₛᵏ
    @test mass_matrix(sst_rand) == Mᵏ
    @test damping_matrix(sst_rand) == Cᵏ
    @test stiffness_matrix(sst_rand) == Kᵏ
    @test strain(sst_rand) == ϵᵏ
    @test stress(sst_rand) == σᵏ

    reset!(sst_rand)
    @test iszero(displacements(sst_rand))
    @test iszero(velocity(sst_rand))
    @test iszero(acceleration(sst_rand))
    @test iszero(Δ_displacements(sst_rand))
    @test iszero(external_forces(sst_rand))
    @test iszero(internal_forces(sst_rand))
    @test iszero(viscous_forces(sst_rand))
    @test iszero(inertial_forces(sst_rand))
    @test iszero(tangent_matrix(sst_rand))
    @test iszero(mass_matrix(sst_rand))
    @test iszero(damping_matrix(sst_rand))
    @test iszero(stiffness_matrix(sst_rand))
    map(strain(sst_rand)) do ϵ
        @test iszero(ϵ)
    end
    map(stress(sst_rand)) do ϵ
        @test iszero(ϵ)
    end
end

@testset "ONSAS.StructuralAnalyses.DynamicAnalyses.DynamicState" begin
    @test displacements(sst_rand) == Uᵏ
    @test velocity(sst_rand) == Udotᵏ
    @test acceleration(sst_rand) == Udotdotᵏ
    @test strain(sst_rand) == ϵᵏ
    @test stress(sst_rand) == σᵏ
end
