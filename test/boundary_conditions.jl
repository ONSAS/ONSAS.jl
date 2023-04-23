###################################
# BoundaryConditions module tests #
###################################
using Test
using ONSAS.BoundaryConditions
using ONSAS.Elements

# Entities 
n₁ = Node(0, 0, 0, dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)], :T => [Dof(25)]]))
n₂ = Node(0, 1, 0, dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)], :T => [Dof(26)]]))
n₃ = Node(0, 0, 1, dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)], :T => [Dof(27)]]))
n₄ = Node(2, 0, 1, dictionary([:u => [Dof(10), Dof(11), Dof(12)], :θ => [Dof(22), Dof(23), Dof(24)], :T => [Dof(28)]]))

t_face = TriangularFace(n₁, n₂, n₃)
tetra = Tetrahedron(n₁, n₂, n₃, n₄)


@testset "ONSAS.BoundaryConditions.FixedDofBoundaryCondition" begin

    # Generic labeled boundary condition
    generic_fixed_dofs = [:u, :θ]
    fixed_components = [1, 3] # Fixes first, second and third component of the corresponding dofs
    generic_bc_label = :fixed_bc_generic
    fixed_bc = FixedDofBoundaryCondition(generic_fixed_dofs, fixed_components, generic_bc_label)

    @test dofs(fixed_bc) == generic_fixed_dofs
    @test components(fixed_bc) == fixed_components
    @test label(fixed_bc) == generic_bc_label

    @test _apply(fixed_bc, n₃) == [Dof(7), Dof(9), Dof(19), Dof(21)]
    @test _apply(fixed_bc, t_face) == [Dof(1), Dof(3), Dof(13), Dof(15), Dof(4),
        Dof(6), Dof(16), Dof(18), Dof(7), Dof(9), Dof(19), Dof(21)]
    @test _apply(fixed_bc, tetra) == [Dof(1), Dof(3), Dof(13), Dof(15), Dof(4), Dof(6), Dof(16),
        Dof(18), Dof(7), Dof(9), Dof(19), Dof(21), Dof(10), Dof(12), Dof(22), Dof(24)]

end

t_to_test = 2.0

@testset "ONSAS.BoundaryConditions.GlobalLoadBoundaryCondition" begin

    # Generic labeled global load boundary condition
    dofs_to_apply_bc = [:u, :θ]
    load_fact_generic(t) = [sin(t), t, t^2]
    generic_values = t -> load_fact_generic(t) .* [1, 1, 1]
    generic_bc_label = :bc_generic
    generic_bc = GlobalLoadBoundaryCondition(
        dofs_to_apply_bc, generic_values, generic_bc_label
    )

    @test dofs(generic_bc) == dofs_to_apply_bc
    @test values(generic_bc) == generic_values
    @test label(generic_bc) == generic_bc_label
    @test generic_bc(t_to_test) == values(generic_bc)(t_to_test)

    # Node force computation 
    loaded_dofs, f_vec = _apply(generic_bc, n₁, t_to_test)
    @test loaded_dofs == [Dof(1), Dof(2), Dof(3), Dof(13), Dof(14), Dof(15)]
    @test f_vec == repeat(generic_bc(t_to_test), 2)

    # Face tension computation 
    loaded_dofs, p_vec = _apply(generic_bc, t_face, t_to_test)
    @test loaded_dofs == vcat(Dof.(1:9), Dof.(13:21))
    @test p_vec == repeat(generic_bc(t_to_test) * area(t_face) / 3, 6)

    # Volume tension computation 
    loaded_dofs, b_vec = _apply(generic_bc, tetra, t_to_test)
    @test loaded_dofs == vcat(Dof.(1:24))

    @test b_vec == repeat(generic_bc(t_to_test) * volume(tetra) / 4, 8)

end

@testset "ONSAS.BoundaryConditions.LocalPressureBoundaryCondition" begin

    # Generic labeled global load boundary condition
    dofs_to_apply_bc = [:u]
    load_fact_generic(t) = t^2
    generic_values = t -> load_fact_generic(t) .* [1]
    generic_bc_label = :bc_generic
    generic_bc = LocalPressureBoundaryCondition(dofs_to_apply_bc, generic_values)

    @test dofs(generic_bc) == dofs_to_apply_bc
    @test values(generic_bc) == generic_values
    @test generic_bc(t_to_test) == values(generic_bc)(t_to_test)

    # Face tension computation 
    loaded_dofs, p_vec = _apply(generic_bc, t_face, t_to_test)
    @test loaded_dofs == vcat(Dof.(1:9))
    @test p_vec == repeat([generic_values(t_to_test)[1] * area(t_face) / 3, 0, 0], 3)
end

