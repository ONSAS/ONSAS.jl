using Test, Dictionaries

using ONSAS.Nodes

@testset "ONSAS.Dof" begin

    # Default dof
    index_θⱼ = 2
    θⱼ_dof = Dof(index_θⱼ)
    @test Nodes.index(θⱼ_dof) == index_θⱼ

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
x_sa2D = Point(xᵢ, 2xᵢ)
x_vec2D = [xᵢ, 2xᵢ]
x_tup2D = (xᵢ, 2xᵢ)
x_test_vec_2D = [x_sa2D, x_vec2D, x_tup2D]
x_test_2D = rand(x_test_vec_2D)
x_sa3D = Point(xᵢ, 2xᵢ, 3xᵢ)
x_vec3D = [xᵢ, 2xᵢ, 3xᵢ]
x_tup3D = (xᵢ, 2xᵢ, 3xᵢ)
x_test_vec_3D = [x_sa3D, x_vec3D, x_tup3D]
x_test_3D = rand(x_test_vec_3D)

@testset "ONSAS.Node" begin
    L = 1

    # Node with empty dofs constructors
    n = Node(L, L, L)
    n = Node((L, L, L))
    n = Node([L, L, L])

    @test coordinates(n) == n.x
    @test isempty(dofs(n))
    @test dimension(n) == length(coordinates(n))
    @test n[1] == coordinates(n)[1]

    # Add dofs
    field = :u
    dofs_indexes = [1, 2, 3, 4]
    apply!(n, field, dofs_indexes)
    @test dofs(n, field) == dofs_indexes

    # Constructor with dofs 
    dofs_dict = dictionary([field => dofs_indexes])
    n_dofs = Node((L, L, L), dofs_dict)
    @test dofs(n_dofs, field) == dofs_indexes
end

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
