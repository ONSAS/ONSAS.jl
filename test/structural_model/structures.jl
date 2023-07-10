using Test
using Dictionaries: dictionary
using ONSAS.Structures
using ONSAS.Squares
using ONSAS.SVKMaterial
using ONSAS.FixedDofBoundaryConditions
using ONSAS.GlobalLoadBoundaryConditions
using ONSAS.Structures
using ONSAS.StructuralEntities
using ONSAS.StructuralMaterials
using ONSAS.StructuralBoundaryConditions
using ONSAS.Nodes
using ONSAS.TriangularFaces
using ONSAS.Trusses
using ONSAS.Tetrahedrons

# Scalar parameters
d = 0.1
L = 2.0
h = 1.0
E = 2e9
ν = 0.3
# Nodes
n₁ = Node(0, 0, 0,
          dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)],
                      :T => [Dof(25)]]))
n₂ = Node(0, 1, 0,
          dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)],
                      :T => [Dof(26)]]))
n₃ = Node(0, 0, 1,
          dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)],
                      :T => [Dof(27)]]))
n₄ = Node(1, 1, 1,
          dictionary([:u => [Dof(10), Dof(11), Dof(12)], :θ => [Dof(22), Dof(23), Dof(24)],
                      :T => [Dof(28)]]))
# Faces
face₁ = TriangularFace(n₁, n₂, n₃)
face₂ = TriangularFace(n₃, n₄, n₃)
# Cross section
s = Square(d)
# Entities
truss₁ = Truss(n₁, n₂, s)
truss₂ = Truss(n₂, n₃, s)
truss₃ = Truss(n₁, n₃, s)
truss₄ = Truss(n₄, n₃, s)

# Materials
steel = SVK(E, ν, "steel")
new_steel = SVK(2E, ν, "new_steel")
aluminum = SVK(E / 3, ν, "aluminium")
mat_dict = dictionary([steel => [truss₁, truss₃], aluminum => [truss₂]])
s_materials = StructuralMaterial(mat_dict)

empty_mat_dict = dictionary([steel=>Vector{AbstractElement}(),aluminum => Vector{AbstractElement}()])
empty_materials = StructuralMaterial(matdict)

@testset "ONSAS.StructuralMaterial" begin
    @test s_materials[truss₁] == steel
    @test s_materials["steel"] == steel
    @test truss₁ ∈ s_materials[steel] && truss₃ ∈ s_materials[steel]
    insert!(s_materials, new_steel)
    @test new_steel ∈ keys(element_materials(s_materials))
    delete!(s_materials, new_steel)
    @test new_steel ∉ keys(element_materials(s_materials))
    material_label_to_be_replaced = "steel"
    replace!(s_materials, new_steel, material_label_to_be_replaced)
    @test new_steel ∈ keys(element_materials(s_materials))
    @test steel ∉ keys(element_materials(s_materials))
    @test s_materials[truss₁] == new_steel
    @test all(map(isempty,element_materials(empty_materials)))
end

# Boundary conditions
Fⱼ = 20.0
Fᵢ = 10.0
dof_dim = 3
bc₁ = FixedDof(:u, collect(1:dof_dim), "fixed_uₓ_uⱼ_uₖ")
bc₂ = FixedDof(:u, [2], "fixed_uⱼ")
bc₃ = GlobalLoad(:u, t -> [0, Fⱼ * t, 0], "load in j")
bc₄ = GlobalLoad(:u, t -> [Fᵢ * sin(t), 0, 0], "load in i")
bc₅ = FixedDof(:T, [1], "fixed_T")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂, n₁]])
face_bc = dictionary([bc₃ => [face₁], bc₅ => [face₁]])
elem_bc = dictionary([bc₄ => [truss₁, truss₂]])

s_boundary_conditions_only_nodes = StructuralBoundaryCondition(; node_bcs=node_bc)
s_boundary_conditions_only_faces = StructuralBoundaryCondition(; face_bcs=face_bc)
s_boundary_conditions_only_elements = StructuralBoundaryCondition(; element_bcs=elem_bc)
s_boundary_conditions = StructuralBoundaryCondition(node_bc, face_bc, elem_bc)

@testset "ONSAS.StructuralBoundaryCondition" begin

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
    @test face₁ ∈ s_boundary_conditions[bc₃] && n₂ ∈ s_boundary_conditions[bc₃] &&
          n₁ ∈ s_boundary_conditions[bc₃]
    @test length(s_boundary_conditions[bc₃]) == 3
    @test bc₂ ∈ s_boundary_conditions[n₂] && bc₃ ∈ s_boundary_conditions[n₂]
    @test bc₄ ∈ s_boundary_conditions[truss₁]

    # Constructor only with node or element boundary onditions(node_bc)
    @test isempty(element_bcs(s_boundary_conditions_only_nodes)) &&
          isempty(face_bcs(s_boundary_conditions_only_nodes))
    @test isempty(element_bcs(s_boundary_conditions_only_faces)) &&
          isempty(node_bcs(s_boundary_conditions_only_faces))

    # Apply boundary conditions
    @test apply(s_boundary_conditions, bc₁) == vcat(Dof.(1:3), Dof.(7:9))
    @test apply(s_boundary_conditions, bc₂) == [Dof(5)]
    @test apply(s_boundary_conditions, bc₅) == Dof.(25:27)

    t_to_test = first(rand(1))
    dofs_bc₃_nodes, f_bc₃_nodes = apply(s_boundary_conditions_only_nodes, bc₃, t_to_test)
    dofs_load_nodes_bc₃ = dictionary(dofs_bc₃_nodes .=> f_bc₃_nodes)
    # dofs__bc₃_nodes_to_test = vcat(dofs(n₁)[dofs(bc₃)...], dofs(n₂)[dofs(bc₃)...])
    dofs_bc₃_faces, f_bc₃_faces = apply(s_boundary_conditions_only_faces, bc₃, t_to_test)
    dofs_load_faces_bc₃ = dictionary(dofs_bc₃_faces .=> f_bc₃_faces)

    dofs_load_bc₃ = mergewith(+, dofs_load_nodes_bc₃, dofs_load_faces_bc₃)
    dofs_bc₃_to_test = collect(keys(dofs_load_bc₃))
    f_bc₃_to_test = collect(values(dofs_load_bc₃))

    dofs_bc₃, f_bc₃ = apply(s_boundary_conditions, bc₃, t_to_test)

    @test dofs_bc₃_to_test == dofs_bc₃
    @test f_bc₃_to_test == f_bc₃

    # Push an entity to a boundary condition
    push!(s_boundary_conditions, bc₃, n₁)
    push!(s_boundary_conditions, bc₃, face₂)
    push!(s_boundary_conditions, bc₄, truss₄)

    @test n₁ ∈ s_boundary_conditions[bc₃] &&
          face₂ ∈ s_boundary_conditions[bc₃] &&
          truss₄ ∈ s_boundary_conditions[bc₄]

    scale_factor = 10
    old_label_to_replace = label(bc₃)
    old_nodes_bc = node_bcs(s_boundary_conditions)[bc₃]
    old_faces_bc = face_bcs(s_boundary_conditions)[bc₃]

    new_bc = GlobalLoad(:u, t -> scale_factor * [0, Fⱼ * t, 0], old_label_to_replace)

    replace!(s_boundary_conditions, new_bc)
    t_to_test = rand()
    @test values(s_boundary_conditions[old_label_to_replace])(t_to_test) ==
          scale_factor * [0, Fⱼ * t_to_test, 0]
    @test node_bcs(s_boundary_conditions)[new_bc] == old_nodes_bc
    @test face_bcs(s_boundary_conditions)[new_bc] == old_faces_bc
end

@testset "ONSAS.StructuralEntity" begin
    sec = Square(1)

    tetra_label = "tetra_label"
    truss_label = "truss_label"
    face_label = "triangle_label"

    velems = [Tetrahedron(tetra_label), Truss(sec, truss_label)]
    vfaces = [TriangularFace(face_label)]

    s_entities = StructuralEntity(velems, vfaces)

    @test elem_types_to_elements(s_entities) == s_entities.elem_types_to_elements
    @test face_types_to_faces(s_entities) == s_entities.face_types_to_faces
    @test elem_types(s_entities) == velems
    @test face_types(s_entities) == vfaces
    @test all([e ∈ all_entities(s_entities) for e in vcat(velems, vfaces)])
    @test s_entities[tetra_label] == velems[1]
    @test s_entities[truss_label] == velems[2]
end

@testset "ONSAS.Structure" begin
    n₁ = Node(0, 0, 0)
    n₂ = Node(0, 1, 0)
    n₃ = Node(0, 0, 1)

    s_mesh = Mesh(; nodes=[n₁, n₂, n₃], elements=[truss₁, truss₂, truss₃])
    set_dofs!(s_mesh, :u, dof_dim)
    s = Structure(s_mesh, s_materials, s_boundary_conditions)

    # Dofs
    @test dofs(s) == dofs(mesh(s))
    @test num_dofs(s) == 3 * num_nodes(s)

    # Nodes
    @test nodes(s) == nodes(mesh(s))
    @test num_nodes(s) == length(nodes(mesh(s)))

    # Faces
    @test faces(s) == faces(mesh(s))
    @test num_faces(s) == length(faces(mesh(s)))

    # Entities
    @test elements(s) == elements(mesh(s))
    @test num_elements(s) == length(elements(mesh(s)))

    # Mesh
    @test mesh(s) == s_mesh

    # Materials
    @test materials(s) == s_materials
    # Replace material
    label_material_to_be_replaced = "new_steel"
    replace!(s, steel, label_material_to_be_replaced)
    @test steel ∈ keys(element_materials(materials(s)))
    @test new_steel ∉ keys(element_materials(materials(s)))

    # Boundary conditions
    @test boundary_conditions(s) == s_boundary_conditions

    @test Dof(4) ∈ free_dofs(s) && Dof(6) ∈ free_dofs(s) && length(free_dofs(s)) == 2
end
