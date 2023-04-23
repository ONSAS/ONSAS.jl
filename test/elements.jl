#########################
# Elements module tests #
#########################
using Test, LinearAlgebra, StaticArrays
using ONSAS.Elements

const RTOL = 1e-3

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
    new_dofs = Dof.(first_dof:(last_dof - 1))
    apply!(node, :u, new_dofs)
    more_new_dofs = Dof.((first_dof + 1):last_dof)
    apply!(node, :u, more_new_dofs)
    new_dofs_node = Dof.(first_dof:last_dof)
    @test length(dofs(node)[:u]) == length(new_dofs_node)
end

@testset "ONSAS.Elements.Node 3D" begin
    node = Node(x_test_3D[1], x_test_3D[2], x_test_3D[3])
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test dimension(node) == length(x_test_3D)
end

@testset "ONSAS.Elements.TriangularFace 3D" begin
    x₁ = [0, 0, 0]
    x₂ = [1, 0, 0]
    x₃ = [0, 1, 0]

    n₁ = Node(x₁, dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)]]))
    n₂ = Node(x₂, dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)]]))
    n₃ = Node(x₃, dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)]]))

    face_label = "my_face"
    f_empty_nodes = TriangularFace(face_label)
    @test label(f_empty_nodes) == Symbol(face_label)
    f₁ = TriangularFace(n₁, n₂, n₃, face_label)
    f₁_no_label = TriangularFace(n₁, n₂, n₃)
    f₁ = create_entity(f_empty_nodes, [n₁, n₂, n₃])

    @test all([n ∈ nodes(f₁) for n in [n₁, n₂, n₃]])
    @test coordinates(f₁) == [coordinates(n₁), coordinates(n₂), coordinates(n₃)]
    @test dimension(f₁) == length(x₁)
    @test all([d ∈ dofs(f₁)[:u] for d in Dof.(1:9)])
    @test all([d ∈ dofs(f₁)[:θ] for d in Dof.(13:21)])
    @test label(f₁) == Symbol(face_label)
    @test area(f₁) == 0.5
    @test normal_direction(f₁) == [0, 0, 1]

    # create entity for gmsh
    empty_entity = TriangularFace(face_label)
    tf = create_entity(empty_entity, [n₁, n₂, n₃])
    @test all([n ∈ nodes(tf) for n in [n₁, n₂, n₃]])
    @test coordinates(tf) == [coordinates(n₁), coordinates(n₂), coordinates(n₃)]
    @test label(empty_entity) == label(tf)
end

E = 1.0
ν = 0.3
my_svk_mat = SVK(; E=E, ν=ν)

@testset "ONSAS.Elements.Truss 1D" begin

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

    t = Truss(n₁, n₂, circle_cross_section, my_label)
    t_no_label = Truss(n₁, n₂, circle_cross_section)
    t_empty_nodes = Truss(circle_cross_section, my_label)
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

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_svk_mat, t, u_global_structure[local_dofs(t)])
    ϵ_rot_ing = (l_def^2 - l_ref^2) / (l_ref * (l_ref + l_def))

    @test ϵ_e[1, 1] ≈ ϵ_rot_ing rtol = RTOL
    @test σ_e[1, 1] ≈ E * ϵ_rot_ing rtol = RTOL
    @test fᵢₙₜ_e[1] ≈ -E * ϵ_rot_ing * A rtol = RTOL
    @test fᵢₙₜ_e[2] ≈ E * ϵ_rot_ing * A rtol = RTOL
    @test Kᵢₙₜ_e ≈ Kᵢₙₜ_e[1, 1] * [1 -1; -1 1] rtol = RTOL
end

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
    l_def = norm(x₂ + u_global_structure[[Dof(1), Dof(3), Dof(5)]] -
                 (x₁ + u_global_structure[[Dof(7), Dof(9), Dof(11)]]))

    A = 1
    square_cross_section = Square(A)
    my_label = "my_3D_truss"
    t = Truss(n₁, n₂, square_cross_section, my_label)
    t_no_label = Truss(n₁, n₂, square_cross_section)
    t_empty_nodes = Truss(square_cross_section, my_label)
    t = create_entity(t_empty_nodes, [n₁, n₂])
    @test label(t_empty_nodes) == Symbol(my_label)

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

n₁ = Node(0, 0, 0, dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)]]))
n₂ = Node(0, 1, 0, dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)]]))
n₃ = Node(0, 0, 1, dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)]]))
n₄ = Node(2, 0, 1,
          dictionary([:u => [Dof(10), Dof(11), Dof(12)], :θ => [Dof(22), Dof(23), Dof(24)]]))

λ = 0.5769
G = 0.3846
my_svk_mat = SVK(λ, G)

tetra_label = "my_tetrahedron"
tetra = Tetrahedron(n₁, n₂, n₃, n₄, tetra_label)

# Global displacements vector of the nodes 
u_global₁_u = [0.1, 0.2, 0.3]
u_global₁_θ = rand(3)
u_global₂_u = [0.4, 0.5, 0.6]
u_global₂_θ = rand(3)
u_global₃_u = [0.7, 0.8, 0.9]
u_global₃_θ = rand(3)
u_global₄_u = [1.0, 1.1, 1.2]
u_global₄_θ = rand(3)

u_global_structure = vcat(u_global₁_u, u_global₂_u, u_global₃_u, u_global₄_u,
                          u_global₁_θ, u_global₂_θ, u_global₃_θ, u_global₄_θ)
n₁ = Node(0.0, 0.0, 0.0,
          dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)]]))
n₂ = Node(0.0, 1.0, 0.0,
          dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)]]))
n₃ = Node(0.0, 0.0, 1.0,
          dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)]]))
n₄ = Node(2.0, 0.0, 1.0,
          dictionary([:u => [Dof(10), Dof(11), Dof(12)], :θ => [Dof(22), Dof(23), Dof(24)]]))

@testset "ONSAS.Elements.Tetrahedron 3D SVK" begin
    tetra_no_label = Tetrahedron(n₁, n₂, n₃, n₄)
    tetra_empty_nodes = Tetrahedron(tetra_label)
    @test label(tetra_empty_nodes) == Symbol(tetra_label)
    tetra = create_entity(tetra_empty_nodes, [n₁, n₂, n₃, n₄])

    @test length(nodes(tetra)) == 4
    @test all([n ∈ nodes(tetra) for n in [n₁, n₂, n₃, n₄]])
    @test all([n ∈ coordinates(tetra) for n in coordinates([n₁, n₂, n₃, n₄])])
    tetra_dofs = dofs(tetra)
    @test all([d ∈ tetra_dofs[:u] for d in Dof.(1:12)])
    @test length(tetra_dofs[:u]) == 12
    @test all([d ∈ tetra_dofs[:θ] for d in Dof.(13:24)])
    @test length(tetra_dofs[:θ]) == 12
    @test local_dof_symbol(tetra) == [:u]
    local_dofs(tetra)
    @test all([d ∈ local_dofs(tetra) for d in Dof.(1:12)])

    @test volume(tetra) == 2 * 1 / 6

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_svk_mat, tetra,
                                               u_global_structure[local_dofs(tetra)])

    # Values from ONSAS.m
    fᵢₙₜ_e_test = [-0.9160, -1.3446, -1.5253, 0.3319, 0.7067, 0.4415,
                   0.3120, 0.5210, 0.9390, 0.2720, 0.1169, 0.1448]

    Kᵢₙₜ_e_test = [           2.1635e+00 7.8458e-01 8.6150e-01 -9.4812e-01 -4.1633e-01 -2.8172e-01 -9.4668e-01 -2.1522e-01 -4.2675e-01 -2.6874e-01 -1.5304e-01 -1.5304e-01
                   7.8458e-01 3.1379e+00 1.5089e+00 -3.7787e-01 -1.6917e+00 -6.0222e-01 -1.8797e-01 -1.2976e+00 -8.6102e-01 -2.1874e-01 -1.4855e-01 -4.5671e-02
                   8.6150e-01 1.5089e+00 3.2917e+00 -3.0095e-01 -7.2401e-01 -1.1596e+00 -3.4181e-01 -7.3923e-01 -1.9835e+00 -2.1874e-01 -4.5671e-02 -1.4855e-01
                   -9.4812e-01 -3.7787e-01 -3.0095e-01 7.0582e-01 2.4326e-01 1.8557e-01 1.4951e-01 3.4454e-02 8.8939e-02 9.2785e-02 1.0016e-01 2.6441e-02
                   -4.1633e-01 -1.6917e+00 -7.2401e-01 2.4326e-01 1.2571e+00 3.0095e-01 2.6441e-02 3.6585e-01 4.0143e-01 1.4663e-01 6.8747e-02 2.1634e-02
                   -2.8172e-01 -6.0222e-01 -1.1596e+00 1.8557e-01 3.0095e-01 8.2120e-01 6.0094e-02 2.8444e-01 2.9374e-01 3.6056e-02 1.6826e-02 4.4710e-02
                   -9.4668e-01 -1.8797e-01 -3.4181e-01 1.4951e-01 2.6441e-02 6.0094e-02 8.8031e-01 1.5204e-01 2.0812e-01 -8.3150e-02 9.4948e-03 7.3595e-02
                   -2.1522e-01 -1.2976e+00 -7.3923e-01 3.4454e-02 3.6585e-01 2.8444e-01 1.5204e-01 1.0165e+00 4.7654e-01 2.8725e-02 -8.4752e-02 -2.1754e-02
                   -4.2675e-01 -8.6102e-01 -1.9835e+00 8.8939e-02 4.0143e-01 2.9374e-01 2.0812e-01 4.7654e-01 1.7697e+00 1.2968e-01 -1.6946e-02 -7.9945e-02
                   -2.6874e-01 -2.1874e-01 -2.1874e-01 9.2785e-02 1.4663e-01 3.6056e-02 -8.3150e-02 2.8725e-02 1.2968e-01 2.5910e-01 4.3388e-02 5.3003e-02
                   -1.5304e-01 -1.4855e-01 -4.5671e-02 1.0016e-01 6.8747e-02 1.6826e-02 9.4948e-03 -8.4752e-02 -1.6946e-02 4.3388e-02 1.6456e-01 4.5791e-02
                   -1.5304e-01 -4.5671e-02 -1.4855e-01 2.6441e-02 2.1634e-02 4.4710e-02 7.3595e-02 -2.1754e-02 -7.9945e-02 5.3003e-02 4.5791e-02 1.8379e-01]

    𝔼_e_test = [        1.3675 0.585 1.02
                0.585 1.87 1.44
                1.02 1.44 3.28]

    σ_e_test = [        -5.9378 -7.8126 -9.5331
                1.6136 1.2953 1.6735
                1.7335 2.5078 4.0564]

    @test fᵢₙₜ_e ≈ fᵢₙₜ_e_test rtol = RTOL
    @test Kᵢₙₜ_e ≈ Kᵢₙₜ_e_test rtol = RTOL
    @test 𝔼_e_test ≈ ϵ_e rtol = RTOL
    # @test σ_e_test ≈ σ_e rtol = RTOL skip = true

    # create entity for gmsh
    empty_tetrahedron = Tetrahedron(tetra_label)
    new_tetra = create_entity(empty_tetrahedron, [n₁, n₂, n₃, n₄])

    # Test weights for interpolation
    # at the nodes should be one
    w₁ = weights(tetra, coordinates(n₁))
    w₄ = weights(tetra, coordinates(n₄))
    @test w₁ ≈ [1.0, 0.0, 0.0, 0.0] rtol = RTOL
    @test w₄ ≈ [0.0, 0.0, 0.0, 1.0] rtol = RTOL

    # The interpolation for a linear scalar field shloud be exact 
    scalar_linear_field(x, y, z) = 10x + 20y + 30z + 40
    sol_at_tetra_nodes = [scalar_linear_field(coordinates(n)...) for n in nodes(tetra)]
    p = [0.5, 0.5, 0.5]
    exact_solution = scalar_linear_field(p...)
    interpolated_solution = dot(sol_at_tetra_nodes, weights(tetra, p))
    @test interpolated_solution ≈ exact_solution rtol = RTOL
end

@testset "ONSAS.Elements.Tetrahedron 3D IsotropicLinearElastic" begin
    my_lin_mat = IsotropicLinearElastic(elasticity_modulus(my_svk_mat),
                                        shear_modulus(my_svk_mat))

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_lin_mat, tetra,
                                               u_global_structure[local_dofs(tetra)])

    # Test internal forces with an HyperElastic material model and zero 𝑢
    equivalent_svk = SVK(lame_parameters(my_lin_mat)...)
    _, Kᵢₙₜ_e_svk, A_, B = internal_forces(equivalent_svk, tetra, zeros(12))

    fᵢₙₜ_e_svk = Kᵢₙₜ_e_svk * u_global_structure[local_dofs(tetra)]

    @test fᵢₙₜ_e_svk ≈ fᵢₙₜ_e rtol = RTOL
    @test Kᵢₙₜ_e_svk ≈ Kᵢₙₜ_e rtol = RTOL
end
