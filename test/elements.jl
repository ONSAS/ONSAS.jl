#########################
# Elements module tests #
#########################
using Test: @testset, @test
using LinearAlgebra: norm
using StaticArrays: SVector
using ONSAS.Elements

@testset "ONSAS.Elements.Dof" begin

    # Default dof
    index_θⱼ = 2
    θⱼ_dof = Dof(index_θⱼ)
    @test index(θⱼ_dof) == index_θⱼ

    u = [2, 3, 4]
    @test u[θⱼ_dof] == u[index_θⱼ]
    new_val = 6
    u[θⱼ_dof] = new_val
    @test u[θⱼ_dof] == new_val

end

# Using StaticArrays, Tuples or Vector
node_eltypes = [Float32, Float64, Int]
node_eltype = rand(node_eltypes)
xᵢ = node_eltype(rand(-100:100))
x_sa2D = SVector(xᵢ, 2xᵢ)
x_vec2D = [xᵢ, 2xᵢ]
x_tup2D = (xᵢ, 2xᵢ)
x_test_vec_2D = [x_sa2D, x_vec2D, x_tup2D]
x_test_2D = rand(x_test_vec_2D)
x_sa3D = SVector(xᵢ, 2xᵢ, 3xᵢ)
x_vec3D = [xᵢ, 2xᵢ, 3xᵢ]
x_tup3D = (xᵢ, 2xᵢ, 3xᵢ)
x_test_vec_3D = [x_sa3D, x_vec3D, x_tup3D]
x_test_3D = rand(x_test_vec_3D)

@testset "ONSAS.Elements.Node 2D" begin

    node = Node(x_test_2D[1], x_test_2D[2])
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test dimension(node) == length(x_vec2D)

    # Dofs
    first_dof = 1
    last_dof = 4
    new_dofs = Dof.(first_dof:last_dof-1)
    add!(node, :u, new_dofs)
    more_new_dofs = Dof.(first_dof+1:last_dof)
    add!(node, :u, more_new_dofs)
    new_dofs_node = Dof.(first_dof:last_dof)
    @test length(dofs(node)[:u]) == length(new_dofs_node)

end

@testset "ONSAS.Elements.Node 3D" begin

    node = Node(x_test_3D[1], x_test_3D[2], x_test_3D[3])
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test dimension(node) == length(x_test_3D)

end


@testset "ONSAS.Elements.TriangularFace" begin

    x₁ = [-1, 0, 0]
    x₂ = [0, 1, 0]
    x₃ = [0, 0, 1]

    n₁ = Node((0, 0, 0), dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)]]))
    n₂ = Node((1, 0, 0), dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)]]))
    n₃ = Node((0, 1, 0), dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)]]))

    face_label = "my_face"
    f₁ = TriangularFace(n₁, n₂, n₃, face_label)
    f₁_no_label = TriangularFace(n₁, n₂, n₃)
    @test all([n ∈ nodes(f₁) for n in [n₁, n₂, n₃]])
    @test coordinates(f₁) == [coordinates(n₁), coordinates(n₂), coordinates(n₃)]
    @test dimension(f₁) == length(x₁)
    @test all([d ∈ dofs(f₁)[:u] for d in Dof.(1:9)])
    @test all([d ∈ dofs(f₁)[:θ] for d in Dof.(13:21)])
    @test label(f₁) == Symbol(face_label)
    @test area(f₁) == 0.5
    @test normal_direction(f₁) == [0, 0, 1]

end

E = 1.0
ν = 0.3
my_mat = SVK(E, ν)

@testset "ONSAS.Elements.Truss 3D" begin

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
    l_def = norm(
        x₂ + u_global_structure[[Dof(1), Dof(3), Dof(5)]] -
        (x₁ + u_global_structure[[Dof(7), Dof(9), Dof(11)]])
    )

    A = 1
    square_corss_section = Square(A)
    my_label = "my_truss"
    t = Truss(n₁, n₂, square_corss_section, my_label)
    t_no_label = Truss(n₁, n₂, square_corss_section)

    @test dimension(t) == dimension(n₁) == dimension(n₂)
    @test n₁ ∈ nodes(t) && n₂ ∈ nodes(t)
    @test x₁ ∈ coordinates(t) && x₂ ∈ coordinates(t)
    @test cross_section(t) == square_corss_section
    truss_dofs = dofs(t)
    @test all([d ∈ truss_dofs[:u] for d in [Dof(1), Dof(3), Dof(5), Dof(7), Dof(9), Dof(11)]])
    @test all([d ∈ truss_dofs[:θ] for d in [Dof(2), Dof(4), Dof(6), Dof(8), Dof(10), Dof(12)]])
    @test local_dof_symbol(t) == [:u]
    local_dofs(t)
    @test all([d ∈ local_dofs(t) for d in [Dof(1), Dof(3), Dof(5), Dof(7), Dof(9), Dof(11)]])
    @test string(label(t)) == my_label

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_mat, t, u_global_structure[local_dofs(t)])
    strain(t, u_global_structure[local_dofs(t)]) == ϵ_e
    stress(t, u_global_structure[local_dofs(t)]) == σ_e
    @test norm(fᵢₙₜ_e) == 0
    @test Kᵢₙₜ_e[1] == E * A / l_def
    @test norm(σ_e) == 0
    @test norm(ϵ_e) == 0

end