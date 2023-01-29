########################
# Mesh interface tests #
########################
using Test: @testset, @test
using StaticArrays: SVector
using ONSAS.Meshes
using ONSAS.Meshes: DEFAULT_NODE_INDEX




@testset "ONSAS.Meshes.Mesh" begin

    # Create a triangular mesh with trusses
    L = 1
    n₁ = Node((0, 0), 1)
    n₂ = Node((L, L), 2)
    n₃ = Node((2L, 0), 3)

    nodes = [n₁, n₂, n₃]

    # Struct to be done in the Elements module
    struct Truss
        n₁::Node
        n₂::Node
    end

    # Struct to be done in the Elements module
    struct Frame
        n₁::Node
        n₂::Node
    end

    ln = [n₂]
    loaded_set = build_node_set("loaded_nodes", ln)
    fn = [n₁, n₃]
    fixed_set = build_node_set("fixed_nodes", fn)

    node_sets_test = Dict{String,Set{Int}}(
        loaded_set[1] => loaded_set[2],
        fixed_set[1] => fixed_set[2],
    )


    elements = Vector{Union{Truss,Frame}}([Truss(n₁, n₂), Truss(n₂, n₃), Truss(n₃, n₁), Frame(n₁, n₃)])
    mesh = Mesh(nodes, elements, node_sets_test)

    @test element_types(mesh) == eltype(elements)
    @test dimension(mesh) == dimension(n₁)
    @test coordinates_eltype(mesh) == coordinates_eltype(n₁)
    @test node_sets(mesh) == node_sets_test






end
