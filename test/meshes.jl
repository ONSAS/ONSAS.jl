########################
# Mesh interface tests #
########################
using Test: @testset, @test
using ONSAS.Meshes
using ONSAS.Elements
using ONSAS.CrossSections: Circle
using ONSAS.Utils




@testset "ONSAS.Meshes.Mesh" begin

    # Mesh considering only nodes
    # Nodes
    L = 1
    n₁ = Node((0, 0))
    n₂ = Node((L, L))
    n₃ = Node((2L, 0))
    mesh_nodes = [n₁, n₂, n₃]
    ## Elements connectivity
    elem₁_nodes = [n₁, n₂]
    elem₂_nodes = [n₂, n₃]
    vec_conec_elems = [elem₁_nodes, elem₂_nodes]
    mesh = Mesh(mesh_nodes, vec_conec_elems)

    @test coordinates_eltype(mesh) == coordinates_eltype(n₁)
    @test dimension(mesh) == 2
    @test length(dofs(mesh)) == 9
    @test isempty(elements(mesh))
    @test element_nodes(mesh) == vec_conec_elems
    @test length(element_sets(mesh)) == 0
    @test nodes(mesh) == mesh_nodes
    @test length(node_sets(mesh)) == 0

    # Add elements
    n₄ = Node((3L, L))
    push!(mesh, n₄)
    @test nodes(mesh)[end] == n₄
    s₁ = Circle(0.1)
    t = Truss(s₁)
    push!(mesh, t)
    @test nodes(mesh)[end] == n₄

    # getting elements with NodeIndex and ElementIndex
    @test mesh[NodeIndex(3)] == n₃
    @test mesh[ElementIndex(1)] == t

end
