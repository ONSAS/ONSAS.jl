#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using ONSAS.BoundaryConditions

@testset "ONSAS.BoundaryConditions.FixedDofBoundaryCondition" begin

    # Generic labeled boundary condition
    generic_fixed_dofs = [:u, :θ]
    fixed_components = [1, 2, 3] # Fixes first, second and third component of the corresponding dofs
    generic_bc_label = :fixed_bc_generic
    fixed_bc = FixedDofBoundaryCondition(generic_fixed_dofs, fixed_components, generic_bc_label)

    @test dofs(fixed_bc) == generic_fixed_dofs
    @test components(fixed_bc) == fixed_components
    @test label(fixed_bc) == generic_bc_label

end



@testset "ONSAS.BoundaryConditions.GlobalLoadBoundaryCondition" begin

    # Generic labeled global load boundary condition
    dofs_to_apply_bc = [:u]
    load_fact_generic(t) = [sin(t), t, t^2]
    generic_values = t -> load_fact_generic(t) .* [1, 1, 1]
    generic_bc_label = :bc_generic
    generic_bc = GlobalLoadBoundaryCondition(
        dofs_to_apply_bc, generic_values, generic_bc_label
    )

    @test dofs(generic_bc) == dofs_to_apply_bc
    @test values(generic_bc) == generic_values
    @test label(generic_bc) == generic_bc_label
    @test generic_bc(2.0) == values(generic_bc)(2.0)

end

@testset "ONSAS.BoundaryConditions.LocalLoadBoundaryCondition" begin

    # Generic labeled global load boundary condition
    dofs_to_apply_bc = [:u]
    load_fact_generic(t) = [sin(t), t, t^2]
    generic_values = t -> load_fact_generic(t) .* [1, 1, 1]
    generic_bc_label = :bc_generic
    generic_bc = LocalLoadBoundaryCondition(dofs_to_apply_bc, generic_values
    )

    @test dofs(generic_bc) == dofs_to_apply_bc
    @test values(generic_bc) == generic_values
    @test generic_bc(2.0) == values(generic_bc)(2.0)

end

