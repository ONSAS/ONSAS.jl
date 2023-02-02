##########################
# Structural model tests #
##########################
using Test: @testset, @test
using ONSAS.Materials, ONSAS.Elements, ONSAS.BoundaryConditions


@testset "ONSAS.StructuralModels" begin

    # Scalar properties
    E = 2.0
    ν = 0.3
    a = 0.01
    Fⱼ = 2e3
    # -------------------------------
    # Materials
    # -------------------------------
    steel = SVK(E, ν)
    aluminum = SVK(E / 3, ν)
    steel_label = "steel"
    aluminum_label = "aluminum"
    materials_dict = Dict(steel_label => steel, aluminum_label => aluminum)
    materials = StructuralMaterials(materials_dict)
    @test label(steel) == Symbol(steel_label)
    @test label(aluminum) == Symbol(aluminum_label)
    # -------------------------------
    # Geometries
    # -------------------------------
    ## Cross section
    s₁ = Square(a)
    s₂ = Square(2a)
    # -------------------------------
    # Elements
    # -------------------------------
    elem₁ = Truss(s₁)
    elem₁_label = "truss_s₁"
    elem₂ = Truss(s₂)
    elem₂_label = "truss_s₂"
    elements_dict = Dict{String,AbstractElement}(
        elem₁_label => elem₁,
        elem₂_label => elem₂
    )
    elements = StructuralElements(elements_dict)
    @test label(elem₁) == Symbol(elem₁_label)
    @test label(elem₂) == Symbol(elem₂_label)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    bc₁ = FixedDisplacementBoundaryCondition()
    bc_label₁ = "fixed"
    bc₂ = FⱼLoadBoundaryCondition(Fⱼ)
    bc_label₂ = "load"
    boundary_conditions = Dict(
        bc_label₁ => bc₁,
        bc_label₂ => bc₂
    )
    bcs = StructuralBoundaryConditions(boundary_conditions)
    @test label(elem₁) == Symbol(elem₁_label)
    @test label(elem₂) == Symbol(elem₂_label)


end