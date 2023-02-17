########################
# Mesh interface tests #
########################
using Test: @testset, @test
using ONSAS.Meshes
using ONSAS.Utils

@testset "ONSAS.Meshes.Mesh" begin

    # Mesh considering only nodes
    # Nodes
    L = 1
    n₁ = Node(0, 0, 0)
    n₂ = Node(L, L, L)
    n₃ = Node(2L, 0, 5L)
    n₄ = Node(2L, 0, 5L)
    vec_nodes = [n₁, n₂, n₃, n₄]

    ## Elements
    d = 0.2
    s₁ = Circle(d)
    t₁ = Truss(n₁, n₂, s₁)
    t₂ = Truss(n₂, n₃, s₁)
    t₃ = Truss(n₃, n₄, s₁)
    vec_elements = [t₁, t₂, t₃]

    # Sets
    set_nodes = dictionary(["set₂" => Set([1, 2, 3]), "set₁" => Set([3, 4])])
    set_elements = dictionary(["set₂" => Set([1, 2]), "set₁" => Set([1, 2])])

    # Constructors
    mesh_with_sets = Mesh(vec_nodes, vec_elements, set_nodes, set_elements)
    mesh = Mesh(vec_nodes, vec_elements, set_nodes)

    @test dimension(mesh) == dimension(n₁)
    @test all(isempty.(dofs(mesh)))
    @test elements(mesh) == vec_elements
    @test nodes(mesh) == vec_nodes
    @test element_sets(mesh_with_sets) == set_elements
    @test node_sets(mesh_with_sets) == set_nodes

    # Add nodes and elements
    new_node₁ = Node(3L, 0, 5L)
    new_node₂ = Node(3L, 0, 5L)
    push!(mesh, [new_node₁, new_node₂])
    @test last(nodes(mesh)) == new_node₂
    new_t₄ = Truss(n₃, n₄, s₁)
    new_t₅ = Truss(n₂, n₁, s₁)
    push!(mesh, [new_t₄, new_t₅])
    @test last(elements(mesh)) == new_t₅
    # Add dofs 
    u_dim = 3
    add!(mesh, :u, u_dim)
    @test num_nodes(mesh) * u_dim == num_dofs(mesh)
    θ_dim = 3
    add!(mesh, :θ, θ_dim)
    @test num_nodes(mesh) * (u_dim + θ_dim) == num_dofs(mesh)

end
