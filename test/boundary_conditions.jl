#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using ONSAS.BoundaryConditions
using ONSAS.BoundaryConditions: DEFAULT_LABEL

@testset "ONSAS.BoundaryConditions.DisplacementBoundaryCondition" begin

    # Generic labeled boundary condition
    generic_fixed_dofs = [:uᵢ, :uⱼ, :θⱼ, :uₖ, :θₖ]
    generic_fixed_values = zeros(length(generic_fixed_dofs))
    generic_bc_label = :bc_generic
    generic_bc = DisplacementBoundaryCondition(generic_fixed_dofs, generic_fixed_values, generic_bc_label)

    # Accessors 
    @test dofs(generic_bc) == generic_fixed_dofs
    @test values(generic_bc) == generic_fixed_values
    @test label(generic_bc) == generic_bc_label

    # Generic no labeled boundary condition
    no_label_generic_bc = DisplacementBoundaryCondition(dofs=generic_fixed_dofs, values=generic_fixed_values)
    @test dofs(no_label_generic_bc) == generic_fixed_dofs
    @test values(no_label_generic_bc) == generic_fixed_values
    @test label(no_label_generic_bc) == DEFAULT_LABEL

    # Fixed boundary condition
    fixed_bc = FixedDisplacementBoundaryCondition()
    @test dofs(fixed_bc) == [:uᵢ, :θᵢ, :uⱼ, :θⱼ, :uₖ, :θₖ]
    @test values(fixed_bc) == zeros(length(dofs(fixed_bc)))
    @test label(fixed_bc) == DEFAULT_LABEL

    # Pinned boundary condition
    label_pinned = :my_pinned_bc
    pinned_bc = PinnedDisplacementBoundaryCondition(label_pinned)
    @test dofs(pinned_bc) == [:uᵢ, :uⱼ, :uₖ]
    @test values(pinned_bc) == zeros(length(dofs(pinned_bc)))
    @test label(pinned_bc) == label_pinned



end