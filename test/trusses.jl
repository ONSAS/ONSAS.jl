using Test, Dictionaries, StaticArrays, LinearAlgebra

using ONSAS.SvkMaterial
using ONSAS.Trusses
using ONSAS.Circles
using ONSAS.Squares
using ONSAS.Nodes
using ONSAS.Entities

E = 1.0
ν = 0.3
my_svk_mat = Svk(; E=E, ν=ν)

const RTOL = 1e-3

@testset "ONSAS.Entities.Truss 1D" begin

    # General case considering a mesh with rotations 
    x₁ = [-1]
    x₂ = [1]
    n₁ = Node(x₁, dictionary([:u => [Dof(1)], :T => [Dof(2)]]))
    n₂ = Node(x₂, dictionary([:u => [Dof(3)], :T => [Dof(4)]]))
    # global displacements 
    u_gobal_₁ = [0.1, 273] # uᵢ, Tᵢ (node 1)
    u_gobal_₂ = [0.25, 273] # uᵢ, Tᵢ (node 2)
    u_global_structure = vcat(u_gobal_₁, u_gobal_₂)
    l_ref = norm(x₂ - x₁)
    l_def = norm(x₂ + u_global_structure[[Dof(3)]] - (x₁ + u_global_structure[[Dof(1)]]))

    A = 1
    d = sqrt(4 * A / pi)
    circle_cross_section = Circle(d)
    my_label = "my_1D_truss"
    strain = RotatedEngineeringStrain

    t = Truss(SVector(n₁, n₂), circle_cross_section, strain, my_label)
    t_no_label = Truss(n₁, n₂, circle_cross_section)
    t_empty_nodes = Truss(circle_cross_section, strain, my_label)
    t = create_entity(t_empty_nodes, [n₁, n₂])

    @test n₁ ∈ nodes(t) && n₂ ∈ nodes(t)
    @test all([n ∈ coordinates(t) for n in coordinates([n₁, n₂])])
    @test cross_section(t) == circle_cross_section
    truss_dofs = dofs(t)
    @test all([d ∈ truss_dofs[:u] for d in [Dof(1), Dof(3)]])
    @test all([d ∈ truss_dofs[:T] for d in [Dof(2), Dof(4)]])
    @test local_dof_symbol(t) == [:u]
    @test local_dofs(t) == [Dof(1), Dof(3)]
    @test string(label(t)) == my_label
    @test strain_model(t) == strain

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_svk_mat, t, u_global_structure[local_dofs(t)])
    ϵ_rot_ing = (l_def^2 - l_ref^2) / (l_ref * (l_ref + l_def))

    @test ϵ_e[1, 1] ≈ ϵ_rot_ing rtol = RTOL
    @test σ_e[1, 1] ≈ E * ϵ_rot_ing * l_def / l_ref rtol = RTOL
    @test fᵢₙₜ_e[1] ≈ -E * ϵ_rot_ing * A rtol = RTOL
    @test fᵢₙₜ_e[2] ≈ E * ϵ_rot_ing * A rtol = RTOL
    @test Kᵢₙₜ_e ≈ Kᵢₙₜ_e[1, 1] * [1 -1; -1 1] rtol = RTOL
end

@testset "ONSAS.Entities.Truss 3D" begin

    # General case considering a mesh with rotations 
    x₁ = [-1, 0, 0]
    x₂ = [1, 0, 0]
    n₁ = Node(x₁, dictionary([:u => [Dof(1), Dof(3), Dof(5)], :θ => [Dof(2), Dof(4), Dof(6)]]))
    n₂ = Node(x₂, dictionary([:u => [Dof(7), Dof(9), Dof(11)], :θ => [Dof(8), Dof(10), Dof(12)]]))
    # global displacements 
    u_gobal_₁ = [0, 0, 0, 0, 0, 0]# uᵢ, θᵢ :uⱼ, θⱼ uₖ, θₖ (node 1)
    u_gobal_₂ = [0, 0, 0, 0, 0, 0]# uᵢ, θᵢ :uⱼ, θⱼ uₖ, θₖ (node 2)
    u_global_structure = vcat(u_gobal_₁, u_gobal_₂)
    l_ref = norm(x₂ - x₁)
    l_def = norm(x₂ + u_global_structure[[Dof(1), Dof(3), Dof(5)]] -
                 (x₁ + u_global_structure[[Dof(7), Dof(9), Dof(11)]]))

    strain_model_t = RotatedEngineeringStrain
    A = 1
    square_cross_section = Square(A)
    my_label = "my_3D_truss"

    # Constructors
    t = Truss(SVector(n₁, n₂), square_cross_section, strain_model_t, my_label)
    t_no_label = Truss(n₁, n₂, square_cross_section, strain_model_t)
    t_empty_nodes = Truss(square_cross_section, strain_model_t, my_label)
    t = create_entity(t_empty_nodes, [n₁, n₂])

    # Accessors
    @test strain_model(t) == strain_model_t
    @test label(t_empty_nodes) == Symbol(my_label)
    @test strain_model(t_empty_nodes) == strain_model_t

    @test n₁ ∈ nodes(t) && n₂ ∈ nodes(t)
    @test all([n ∈ coordinates(t) for n in coordinates([n₁, n₂])])
    @test cross_section(t) == square_cross_section
    truss_dofs = dofs(t)
    @test all([d ∈ truss_dofs[:u] for d in [Dof(1), Dof(3), Dof(5), Dof(7), Dof(9), Dof(11)]])
    @test all([d ∈ truss_dofs[:θ] for d in [Dof(2), Dof(4), Dof(6), Dof(8), Dof(10), Dof(12)]])
    @test local_dof_symbol(t) == [:u]
    @test local_dofs(t) == [Dof(1), Dof(3), Dof(5), Dof(7), Dof(9), Dof(11)]
    @test string(label(t)) == my_label

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_svk_mat, t, u_global_structure[local_dofs(t)])
    strain(t, u_global_structure[local_dofs(t)]) == ϵ_e
    stress(t, u_global_structure[local_dofs(t)]) == σ_e
    @test norm(fᵢₙₜ_e) == 0
    @test Kᵢₙₜ_e[1] == E * A / l_def
    @test norm(σ_e) == 0
    @test norm(ϵ_e) == 0
end
