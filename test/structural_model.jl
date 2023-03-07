##########################
# Structural model tests #
##########################
using Test: @testset, @test
using ONSAS.StructuralModel

# Scalar parameters
d = 0.1
L = 2.0
h = 1.0
E = 2e9
ν = 0.3
# Nodes
n₁ = Node(0, 0, 0, dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)], :T => [Dof(25)]]))
n₂ = Node(0, 1, 0, dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)], :T => [Dof(26)]]))
n₃ = Node(0, 0, 1, dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)], :T => [Dof(27)]]))
# Faces 
face₁ = TriangularFace(n₁, n₂, n₃)
# Cross section
s = Square(d)
# Elements
truss₁ = Truss(n₁, n₂, s)
truss₂ = Truss(n₂, n₃, s)
truss₃ = Truss(n₁, n₃, s)
# Mesh
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
bc₅ = FixedDofBoundaryCondition([:T], [1], "fixed_T")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂, n₁]])
face_bc = dictionary([bc₃ => [face₁], bc₅ => [face₁]])
elem_bc = dictionary([bc₄ => [truss₁, truss₂]])

s_boundary_conditions_only_nodes = StructuralBoundaryConditions(node_bcs=node_bc)
s_boundary_conditions_only_faces = StructuralBoundaryConditions(face_bcs=face_bc)
s_boundary_conditions_only_elements = StructuralBoundaryConditions(element_bcs=elem_bc)
s_boundary_conditions = StructuralBoundaryConditions(node_bc, face_bc, elem_bc)

@testset "ONSAS.StructuralModel.StructuralMaterials" begin

    @test s_materials[truss₁] == steel
    @test s_materials["steel"] == steel
    @test truss₁ ∈ s_materials[steel] && truss₃ ∈ s_materials[steel]

end

@testset "ONSAS.StructuralModel.StructuralBoundaryConditions" begin

    # Access and filter boundary conditions
    @test node_bcs(s_boundary_conditions) == node_bc
    @test face_bcs(s_boundary_conditions) == face_bc
    @test element_bcs(s_boundary_conditions) == elem_bc
    @test all(bc ∈ all_bcs(s_boundary_conditions) for bc in [bc₁, bc₂, bc₃, bc₄, bc₅])
    @test length(all_bcs(s_boundary_conditions)) == 5

    @test length(load_bcs(s_boundary_conditions)) == 2
    @test bc₃ ∈ load_bcs(s_boundary_conditions) && bc₄ ∈ load_bcs(s_boundary_conditions)
    @test length(fixed_dof_bcs(s_boundary_conditions)) == 3
    @test all(bc ∈ fixed_dof_bcs(s_boundary_conditions) for bc in [bc₁, bc₂, bc₅])

    @test s_boundary_conditions["fixed_uⱼ"] == bc₂
    @test truss₁ ∈ s_boundary_conditions[bc₄]
    @test face₁ ∈ s_boundary_conditions[bc₃] && n₂ ∈ s_boundary_conditions[bc₃] && n₁ ∈ s_boundary_conditions[bc₃]
    @test length(s_boundary_conditions[bc₃]) == 3
    @test bc₂ ∈ s_boundary_conditions[n₂] && bc₃ ∈ s_boundary_conditions[n₂]
    @test bc₄ ∈ s_boundary_conditions[truss₁]

    # Constructor only with node or element boundary onditions(node_bc)
    @test isempty(element_bcs(s_boundary_conditions_only_nodes)) && isempty(face_bcs(s_boundary_conditions_only_nodes))
    @test isempty(element_bcs(s_boundary_conditions_only_faces)) && isempty(node_bcs(s_boundary_conditions_only_faces))

    # Apply boundary conditions
    @test _apply(s_boundary_conditions, bc₁) == vcat(Dof.(1:3), Dof.(7:9))
    @test _apply(s_boundary_conditions, bc₂) == [Dof(5)]
    @test _apply(s_boundary_conditions, bc₅) == Dof.(25:27)

    t_to_test = first(rand(1))
    dofs_bc₃_nodes, f_bc₃_nodes = _apply(s_boundary_conditions_only_nodes, bc₃, t_to_test)
    dofs_load_nodes_bc₃ = dictionary(dofs_bc₃_nodes .=> f_bc₃_nodes)
    # dofs__bc₃_nodes_to_test = vcat(dofs(n₁)[dofs(bc₃)...], dofs(n₂)[dofs(bc₃)...])
    dofs_bc₃_faces, f_bc₃_faces = _apply(s_boundary_conditions_only_faces, bc₃, t_to_test)
    dofs_load_faces_bc₃ = dictionary(dofs_bc₃_faces .=> f_bc₃_faces)

    dofs_load_bc₃ = mergewith(+, dofs_load_nodes_bc₃, dofs_load_faces_bc₃)
    dofs_bc₃_to_test = collect(keys(dofs_load_bc₃))
    f_bc₃_to_test = collect(values(dofs_load_bc₃))

    dofs_bc₃, f_bc₃ = _apply(s_boundary_conditions, bc₃, t_to_test)

    @test dofs_bc₃_to_test == dofs_bc₃
    @test f_bc₃_to_test == f_bc₃

end

@testset "ONSAS.StructuralModel.StructuralEntities" begin

    vec_elems = [Tetrahedron("tetra_label"), Truss("truss_label")]
    vec_faces = [TriangularFace("triangle_label")]

    s_entities = StructuralEntities(vec_elems, vec_faces)


    @test elem_types_to_elements(s_entities)


end


@testset "ONSAS.StructuralModel.Structure" begin

    n₁ = Node(0, 0, 0)
    n₂ = Node(0, 1, 0)
    n₃ = Node(0, 0, 1)

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

    # Boundary conditions
    @test boundary_conditions(s) == s_boundary_conditions
    @test displacement_bcs(s) == displacement_bcs(s_boundary_conditions)
    @test load_bcs(s) == load_bcs(s_boundary_conditions)

    @test Dof(4) ∈ free_dofs(s) && Dof(6) ∈ free_dofs(s) && length(free_dofs(s)) == 2
end