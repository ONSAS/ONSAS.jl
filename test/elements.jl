########################
# Mesh interface tests #
########################
using Test: @testset, @test
using LinearAlgebra: norm
using StaticArrays: SVector
using ONSAS.Materials: SVK
using ONSAS.Elements
using ONSAS.Elements: _dim_to_nodal_dofs

@testset "ONSAS.Elements.Dof" begin

    # Default dof
    index_θⱼ = 4
    θⱼ_dof = Dof(index_θⱼ)
    @test index(θⱼ_dof) == index_θⱼ
    new_index = index_θⱼ + 1
    setindex!(θⱼ_dof, new_index)
    @test index(θⱼ_dof) == new_index

end

@testset "ONSAS.Elements.Node 2D" begin

    # Using StaticArrays, Tuples or Vector
    node_eltypes = [Float32, Float64, Int]
    node_eltype = rand(node_eltypes)
    xᵢ = node_eltype(rand(-100:100))
    x_sa = SVector(xᵢ, 2xᵢ)
    x_vec = [xᵢ, 2xᵢ]
    x_tup = (xᵢ, 2xᵢ)
    x_test_vec = [x_sa, x_vec, x_tup]
    x_test = rand([x_sa, x_vec, x_tup])

    node = Node(x_test[1], x_test[2])
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test dimension(node) == length(x_test)

    # Dofs
    @test all(vcat(dofs(node), dofs(node))[i] == d for (i, d) in enumerate(dofs([node, node])))
    # Set new index and test dofs indexes
    new_index = 2
    new_global_dof_indexes = 4:6
    setindex!(node, new_index)
    @test all([index(dof) == new_global_dof_indexes[i] for (i, dof) in enumerate(dofs(node))])

end

@testset "ONSAS.Elements.Node 3D" begin

    # Using StaticArrays, Tuples or Vector
    node_eltypes = [Float32, Float64, Int]
    node_eltype = rand(node_eltypes)
    xᵢ = node_eltype(rand(-100:100))
    x_sa = SVector(xᵢ, 2xᵢ, 3xᵢ)
    x_vec = [xᵢ, 2xᵢ, 3xᵢ]
    x_tup = (xᵢ, 2xᵢ, 3xᵢ)
    x_test_vec = [x_sa, x_vec, x_tup]
    x_test = rand([x_sa, x_vec, x_tup])

    node = Node(x_test[1], x_test[2], x_test[3])
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test dimension(node) == length(x_test)

    # Dofs
    @test all(vcat(dofs(node), dofs(node))[i] == d for (i, d) in enumerate(dofs([node, node])))
    # Set new index and test dofs indexes
    new_index = 2
    new_global_dof_indexes = 7:12
    setindex!(node, new_index)
    @test all([index(dof) == new_global_dof_indexes[i] for (i, dof) in enumerate(dofs(node))])

end

@testset "ONSAS.Elements.Truss 3D" begin

    E = 1.0
    ν = 0.3
    my_mat = SVK(E, ν)

    x₁ = [-1, 0, 0]
    x₂ = [1, 0, 0]
    l₀ = norm(x₂ - x₁)
    u₁ = [0, 0, 0, 0, 0, 0]
    u₂ = [0, 0, 0, 0, 0, 0]
    l₁ = norm(x₂ + u₂[1:2:end] - (x₁ + u₁[1:2:end]))

    n₁ = Node(x₁)
    n₂ = Node(x₂)

    t_nodes = [n₁, n₂]
    A = 1
    sq = Square(A)
    my_label = "my_truss"
    t = Truss(t_nodes, sq, my_label)

    @test nodes(t) == t_nodes
    @test length(coordinates(t)) == 6
    @test last(index.(dofs(t))) == 12
    @test string(label(t)) == my_label

    # @test dofs(t) = 

    u_e = vcat(u₁, u₂)

    fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(my_mat, t, u_e)
    @test norm(fᵢₙₜ_e) == 0
    @test Kᵢₙₜ_e[1] == E * A / l₁
    @test norm(σ_e) == 0
    @test norm(ϵ_e) == 0

end