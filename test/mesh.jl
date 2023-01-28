########################
# Mesh interface tests #
########################
using Test: @testset, @test
using StaticArrays: SVector
using ONSAS.Mesh
using ONSAS.Mesh: DEFAULT_NODE_INDEX


@testset "ONSAS.Mesh.Node" begin

    # Using StaticArrays, Tuples or Vector
    node_eltypes = [Float32, Float64, Int]
    node_eltype = rand(node_eltypes)
    xᵢ = rand(node_eltype)
    x_sa = SVector(xᵢ, 2xᵢ, 3xᵢ)
    x_vec = [xᵢ, 2xᵢ, 3xᵢ]
    x_tup = (xᵢ, 2xᵢ, 3xᵢ)
    x_test_vec = [x_sa, x_vec, x_tup]
    x_test = rand([x_sa, x_vec, x_tup])

    node = Node(x_test)
    @test all([node[i] == xᵢ for (i, xᵢ) in enumerate(coordinates(node))])
    @test coordinates_eltype(node) == node_eltype
    @test index(node) == DEFAULT_NODE_INDEX
    @test all([node[i] == x for (i, x) in enumerate(x_test)])
    @test dimension(node) == length(x_test)


end
