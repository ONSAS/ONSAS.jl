##########################
# Structural model tests #
##########################
using Test: @testset, @test
using ONSAS


@testset "ONSAS.StructuralModels" begin

    # Scalar properties
    E = 2.0
    ν = 0.3
    a = 0.01
    Fⱼ = 2e3
    L = 2 # Length in m 
    d = L * cos(deg2rad(65))   # vertical distance in m
    h = L * sin(deg2rad(65))
    # -------------------------------
    # Create mesh
    # -------------------------------
    ## Nodes
    n₁ = Node((0.0, 0.0, 0.0))
    n₂ = Node((d, h, 0.0))
    n₃ = Node((2d, 0.0, 0.0))
    vec_nodes = [n₁, n₂, n₃]
    ## Elements connectivity
    elem₁_nodes = [n₁, n₂]
    elem₂_nodes = [n₂, n₃]
    vec_conec_elems = [elem₁_nodes, elem₂_nodes]
    s_mesh = Mesh(vec_nodes, vec_conec_elems)
    # -------------------------------
    # Materials
    # -------------------------------
    steel = SVK(E, ν)
    aluminum = SVK(E / 3, ν)
    steel_label = "steel"
    aluminum_label = "aluminum"
    materials_dict = Dict(steel_label => steel, aluminum_label => aluminum)
    s_materials = StructuralMaterials(materials_dict)
    # Sets
    m_sets = Dict(
        "steel" => Set([ElementIndex(1)]),
        "aluminum" => Set([ElementIndex(2)]),
    )
    add_set!(s_materials, m_sets)
    # 
    @test label(steel) == Symbol(steel_label)
    @test label(aluminum) == Symbol(aluminum_label)
    @test s_materials["steel"] == steel
    @test steel ∈ materials(s_materials) && aluminum ∈ materials(s_materials)
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
    s_elements = StructuralElements(elements_dict)
    # Sets
    e_sets = Dict(
        "truss_s₁" => Set([ElementIndex(1)]),
        "truss_s₂" => Set([ElementIndex(2)]),
    )
    add_set!(s_elements, e_sets)
    #
    @test label(elem₁) == Symbol(elem₁_label)
    @test label(elem₂) == Symbol(elem₂_label)
    @test s_elements[elem₁_label] == elem₁
    @test elem₁ ∈ elements(s_elements) && elem₂ ∈ elements(s_elements)
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
    s_bcs = StructuralBoundaryConditions(boundary_conditions)
    @test label(elem₁) == Symbol(elem₁_label)
    @test label(elem₂) == Symbol(elem₂_label)
    @test s_bcs[bc_label₁] == bc₁
    @test s_bcs[bc_label₂] == bc₂
    @test bc₁ ∈ disp_bcs(s_bcs) && bc₂ ∈ load_bcs(s_bcs)
    # -------------------------------
    # Create Structure
    # -------------------------------
    s = Structure(s_mesh, s_materials, s_elements, s_bcs)

    @test mesh(s) == s_mesh
    @test s[ElementIndex(1)] == first(elements(mesh(s)))
    @test s[NodeIndex(2)] == n₂
end