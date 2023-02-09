########################
# Mesh interface tests #
########################
using Test: @testset, @test
using StaticArrays: SVector
using ONSAS.Elements
using ONSAS.Elements: _DEFAULT_INDEX_INT, _dim_to_nodal_dofs
using ONSAS.BoundaryConditions: FixedDisplacementBoundaryCondition

@testset "ONSAS.Elements.Dof" begin

    # Default dof
    index_θⱼ = 4
    θⱼ_dof = Dof(index_θⱼ)
    @test index(θⱼ_dof) == index_θⱼ
    new_index = index_θⱼ + 1
    setindex!(θⱼ_dof, new_index)
    @test index(θⱼ_dof) == new_index

end

@testset "ONSAS.Elements.Node" begin

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
    @test index(node) == _DEFAULT_INDEX_INT
    @test all([node[i] == x for (i, x) in enumerate(x_test)])
    @test dimension(node) == length(x_test)

    # Dofs
    @test all(vcat(dofs(node), dofs(node))[i] == d for (i, d) in enumerate(dofs([node, node])))
    # Set new index and test dofs indexes
    new_index = 2
    new_global_dof_indexes = 7:12
    setindex!(node, new_index)
    @test index(node) == new_index
    @test all([index(dof)[] == new_global_dof_indexes[i] for (i, dof) in enumerate(dofs(node))])


    # Boundary conds
    @test length(boundary_conditions(node)) == 0
    fixed_bc = FixedDisplacementBoundaryCondition(:my_fixed_bc)
    pinned_bc = FixedDisplacementBoundaryCondition(:my_pinned_bc)
    push!(node, [fixed_bc, pinned_bc])
    @test length(boundary_conditions(node)) == 2

end

@testset "ONSAS.Elements.Truss" begin



end