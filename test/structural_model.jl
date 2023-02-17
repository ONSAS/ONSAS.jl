##########################
# Structural model tests #
##########################
using Test: @testset, @test
using ONSAS.StructuralModel
using ONSAS.BoundaryConditions


# Scalar parameters
d = 0.1
L = 2.0
h = 1.0
E = 2e9
ν = 0.3

# Nodes
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(L, h, 0.0)
n₃ = Node(2L, 0.0, 0.0)

# Cross section
s = Square(d)

# Elements
truss₁ = Truss(n₁, n₂, s)
truss₂ = Truss(n₂, n₃, s)
truss₃ = Truss(n₁, n₃, s)

# Materials
steel = SVK(E, ν, "steel")
aluminum = SVK(E / 3, ν, "aluminium")
mat_dict = dictionary([steel => [truss₁, truss₃], aluminum => [truss₂]])
s_materials = StructuralMaterials(mat_dict)

# Boundary conditions
Fⱼ = 20.0
Fᵢ = 10.0
dof_dim = 3
bc₁ = FixedDofBoundaryCondition([:u], collect(1:dof_dim), "fixed_uₓ_uⱼ_uₖ")
bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed_uⱼ")
bc₃ = GlobalLoadBoundaryCondition([:u], t -> [0, Fⱼ * t, 0], "load in j")
bc₄ = GlobalLoadBoundaryCondition([:u], t -> [Fᵢ * sin(t), 0, 0], "load in i")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂]])
elem_bc = dictionary([bc₄ => [truss₁, truss₂]])
s_boundary_conditions = StructuralBoundaryConditions(node_bc, elem_bc)

@testset "ONSAS.StructuralModel.StructuralMaterials" begin

    @test s_materials[truss₁] == steel
    @test s_materials["steel"] == steel
    @test truss₁ ∈ s_materials[steel] && truss₃ ∈ s_materials[steel]

end


@testset "ONSAS.StructuralModel.StructuralMaterials" begin


    @test node_bcs(s_boundary_conditions) == node_bc
    @test element_bcs(s_boundary_conditions) == elem_bc

    @test length(load_bcs(s_boundary_conditions)) == 2
    @test bc₃ ∈ load_bcs(s_boundary_conditions) && bc₄ ∈ load_bcs(s_boundary_conditions)
    @test length(fixed_dof_bcs(s_boundary_conditions)) == 2
    @test bc₁ ∈ fixed_dof_bcs(s_boundary_conditions) && bc₂ ∈ fixed_dof_bcs(s_boundary_conditions)

    @test s_boundary_conditions["fixed_uⱼ"] == bc₂
    @test truss₁ ∈ s_boundary_conditions[bc₃] && n₃ ∈ s_boundary_conditions[bc₃]
    @test bc₂ ∈ s_boundary_conditions[n₂] && bc₂ ∈ s_boundary_conditions[n₂] && bc₃ ∈ s_boundary_conditions[n₃]
    @test bc₄ ∈ s_boundary_conditions[truss₁]


    # Constructor only with node or element boundary conditions
    s_boundary_conditions_nodes = StructuralBoundaryConditions(node_bc)
    @test isempty(element_bcs(s_boundary_conditions_nodes))
    Main.@infiltrate

    s_boundary_conditions_element = StructuralBoundaryConditions(elem_bc)

end

@testset "ONSAS.StructuralModel.Structure" begin

    s_mesh = Mesh([n₁, n₂, n₃], [truss₁, truss₂, truss₃])
    add!(s_mesh, :u, dof_dim)
    s = Structure(s_mesh, s_materials, s_boundary_conditions)

    # Dofs
    @test dofs(s) == dofs(mesh(s))
    @test num_dofs(s) == 3 * num_nodes(s)

    # Nodes
    @test nodes(s) == nodes(mesh(s))
    @test num_nodes(s) == length(nodes(mesh(s)))

    # Elements
    @test elements(s) == elements(mesh(s))
    @test num_elements(s) == length(elements(mesh(s)))

    # Mesh
    @test mesh(s) == s_mesh
    @test all([Dof(i) ∈ free_dofs(s) for i in 1:num_dofs(s)])

    # Boundary conditions
    @test boundary_conditions(s) == s_boundary_conditions
    @test displacement_bcs(s) == displacement_bcs(s_boundary_conditions)
    @test load_bcs(s) == load_bcs(s_boundary_conditions)

end