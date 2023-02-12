#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using ONSAS.BoundaryConditions

@testset "ONSAS.BoundaryConditions.DisplacementBoundaryCondition" begin

    # Generic labeled boundary condition
    generic_fixed_dofs = [1, 3, 4, 6]
    generic_fixed_values = zeros(length(generic_fixed_dofs))
    generic_bc_label = :bc_generic
    generic_bc = DisplacementBoundaryCondition(generic_fixed_dofs, generic_fixed_values, generic_bc_label)

    @test dofs(generic_bc) == generic_fixed_dofs
    @test values(generic_bc) == generic_fixed_values
    @test label(generic_bc) == generic_bc_label

    # Generic no labeled boundary condition
    no_label_generic_bc = DisplacementBoundaryCondition(dofs=generic_fixed_dofs, values=generic_fixed_values)
    @test dofs(no_label_generic_bc) == generic_fixed_dofs
    @test values(no_label_generic_bc) == generic_fixed_values
    @test label(no_label_generic_bc) == :no_labelled_bc

    # Fixed boundary condition
    fixed_bc = FixedDisplacementBoundaryCondition()
    @test dofs(fixed_bc) == [1, 2, 3, 4, 5, 6]#[:uᵢ, :θᵢ, :uⱼ, :θⱼ, :uₖ, :θₖ]
    t = 0
    @test values(fixed_bc)(t) == zeros(length(dofs(fixed_bc)))
    @test label(fixed_bc) == :no_labelled_bc

    # Pinned boundary condition
    label_pinned = :my_pinned_bc
    pinned_bc = PinnedDisplacementBoundaryCondition(label_pinned)
    @test dofs(pinned_bc) == [1, 3, 5]#[:uᵢ, :uⱼ, :uₖ]
    @test values(pinned_bc)(t) == zeros(length(dofs(pinned_bc)))
    @test label(pinned_bc) == label_pinned

end


@testset "ONSAS.BoundaryConditions.LoadBoundaryCondition" begin

    ###############
    # Constructors
    ###############

    # Generic labeled global load boundary condition
    generic_load_dofs = [:Mᵢ, :Fⱼ, :Mₖ]
    generic_values = rand(3)
    generic_bc_label = :bc_generic
    load_fact_generic(t) = [sin(t), t, t^2]
    generic_bc = GlobalLoadBoundaryCondition(
        generic_load_dofs, generic_values, load_fact_generic, generic_bc_label
    )


    @test dofs(generic_bc) == generic_load_dofs
    @test values(generic_bc) == generic_values
    @test load_factor_function(generic_bc) == load_fact_generic
    @test label(generic_bc) == generic_bc_label


    # Generic labeled global load boundary condition without load factor function defined
    generic_bc = GlobalLoadBoundaryCondition(dofs=generic_load_dofs, values=generic_values)
    @test load_factor_function(generic_bc) == t -> t
    @test label(generic_bc) == :no_labelled_bc

    # Load boundary condition: Moment along `x` axis 
    Mᵢ_val = rand(Float32, 1)
    Mᵢ_load_factor_func = (t) -> sin(t)
    Mᵢ_bc = MᵢLoadBoundaryCondition(Mᵢ_val..., Mᵢ_load_factor_func)
    @test dofs(Mᵢ_bc) == [2]
    @test values(Mᵢ_bc) == Mᵢ_val
    @test label(Mᵢ_bc) == :no_labelled_bc
    @test load_factor_function(Mᵢ_bc) == Mᵢ_load_factor_func

    # Load boundary condition: Moment along `y` axis 
    Mⱼ_val = rand(Int, 1)
    Mⱼ_bc = MⱼLoadBoundaryCondition(Mⱼ_val...)
    @test dofs(Mⱼ_bc) == [4]
    @test values(Mⱼ_bc) == Mⱼ_val
    @test label(Mⱼ_bc) == :no_labelled_bc
    @test load_factor_function(Mⱼ_bc) == t -> t

    # Load boundary condition: Moment along `z` axis 
    Mₖ_val = rand(Int, 1)
    Mₖ_label = :my_moment_bc
    Mₖ_bc = MₖLoadBoundaryCondition(Mₖ_val..., label=Mₖ_label)
    @test dofs(Mₖ_bc) == [6]
    @test values(Mₖ_bc) == Mₖ_val
    @test label(Mₖ_bc) == Mₖ_label
    @test load_factor_function(Mₖ_bc) == t -> t

    # Load boundary condition: force along `x` axis 
    Fᵢ_val = rand(Float32, 1)
    Fᵢ_load_factor_func = (t) -> sin(t)
    Fᵢ_bc = FᵢLoadBoundaryCondition(Fᵢ_val..., Fᵢ_load_factor_func)
    @test dofs(Fᵢ_bc) == [1]
    @test values(Fᵢ_bc) == Fᵢ_val
    @test label(Fᵢ_bc) == :no_labelled_bc
    @test load_factor_function(Fᵢ_bc) == Fᵢ_load_factor_func

    # Load boundary condition: force along `y` axis 
    Fⱼ_val = rand(Int, 1)
    Fⱼ_bc = FⱼLoadBoundaryCondition(Fⱼ_val...)
    @test dofs(Fⱼ_bc) == [3]
    @test values(Fⱼ_bc) == Fⱼ_val
    @test label(Fⱼ_bc) == :no_labelled_bc
    @test load_factor_function(Fⱼ_bc) == t -> t

    # Load boundary condition: force along `z` axis 
    Fₖ_val = rand(Int, 1)
    Fₖ_label = :Fk_force_bc
    Fₖ_bc = FₖLoadBoundaryCondition(Fₖ_val..., label=Fₖ_label)
    @test dofs(Fₖ_bc) == [5]#[:Fₖ]
    @test values(Fₖ_bc) == Fₖ_val
    @test label(Fₖ_bc) == Fₖ_label
    @test load_factor_function(Fₖ_bc) == t -> t

end