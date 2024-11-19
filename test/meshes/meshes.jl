using Test, Suppressor, Dictionaries
using Test
using Suppressor
using Dictionaries

using ONSAS.Meshes
using ONSAS.StructuralEntities
using ONSAS.Gmsh
using ONSAS.Entities
using ONSAS.Tetrahedrons
using ONSAS.TriangularFaces
using ONSAS.Nodes
using ONSAS.Circles
using ONSAS.Trusses
using ONSAS.TriangularFaces
using ONSAS.Tetrahedrons

const RTOL = 1e-5

@testset "ONSAS.Meshes.Mesh Node + Trusses" begin

    # Mesh considering only nodes
    # Nodes
    L = 1
    n₁ = Node(0, 0, 0)
    n₂ = Node(L, L, L)
    n₃ = Node(2L, 0, 5L)
    n₄ = Node(3L, 0, 5L)
    vec_nodes = [n₁, n₂, n₃, n₄]

    ## Entities
    d = 0.2
    s₁ = Circle(d)
    t₁ = Truss(n₁, n₂, s₁)
    t₂ = Truss(n₂, n₃, s₁)
    t₃ = Truss(n₃, n₄, s₁)
    vec_elements = [t₁, t₂, t₃]

    # Constructors
    mesh = Mesh(; nodes = vec_nodes, elements = vec_elements)

    @test dimension(mesh) == dimension(n₁)
    @test all(isempty.(dofs(mesh)))
    @test node(mesh, 1) == n₁
    @test nodes(mesh) == vec_nodes
    @test elements(mesh) == vec_elements
    @test element(mesh, 1) == t₁

    # Standard form
    nodes_matrix = Matrix{Float64}(undef, dimension(mesh), num_nodes(mesh))
    for (i, n) in enumerate(nodes(mesh))
        nodes_matrix[:, i] = coordinates(n)
    end
    @test node_matrix(mesh) == nodes_matrix
    @test Meshes.connectivity(mesh) == [[1, 2], [2, 3], [3, 4]]

    # Add new nodes and elements.
    new_node₁ = Node(3L, 0, 5L)
    new_node₂ = Node(3L, 0, 5L)
    append!(nodes(mesh), [new_node₁, new_node₂])
    @test last(nodes(mesh)) == new_node₂
    new_t₄ = Truss(n₃, new_node₁, s₁)
    new_t₅ = Truss(n₂, new_node₂, s₁)
    append!(elements(mesh), [new_t₄, new_t₅])
    @test last(elements(mesh)) == new_t₅

    # Add dofs
    u_dim = 3
    set_dofs!(mesh, :u, u_dim)
    @test num_nodes(mesh) * u_dim == num_dofs(mesh)
    θ_dim = 3
    set_dofs!(mesh, :θ, θ_dim)
    @test num_nodes(mesh) * (u_dim + θ_dim) == num_dofs(mesh)
end

uniaxial_mesh_path = joinpath(@__DIR__, "..", "..", "examples", "uniaxial_extension",
    "uniaxial_mesh.jl")
include(uniaxial_mesh_path)

@testset "ONSAS.Meshes.GMSH.MshFile " begin
    Lᵢ = 2.0  # Dimension in x of the box in m
    Lⱼ = 1.0  # Dimension in y of the box in m
    Lₖ = 1.0  # Dimension in z of the box in m

    # Labels
    mat_label = "svkHyper"
    faces_label = "triangle"
    elems_label = "tetrahedron"
    entities_labels = [faces_label, elems_label]
    bc₁_label = "fixed-ux"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    # Load
    bc₄_label = "tension"
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
    # labels
    labels = [mat_label, entities_labels, bc_labels]
    file_name = "uniaxial_test_mesh"
    # Create, load  and delete the mesh
    # Refinement mesh factor
    ms = 0.5
    dir = @__DIR__
    @capture_out begin
        file_name = create_uniaxial_mesh(Lᵢ, Lⱼ, Lₖ, labels, file_name, ms, dir)
    end
    msh_path = joinpath(dir, file_name)
    msh_file = MshFile(msh_path)
    rm(file_name; force = true)

    @test nodes(msh_file) == msh_file.vec_nodes
    @test length(physical_index(msh_file)) == length(connectivity(msh_file))
    @test dimension(msh_file) == dimension(first(nodes(msh_file)))
    @test material_label(msh_file) == ["", "", "", "", "svkHyper"]
    @test entity_label(msh_file) ==
          [faces_label, faces_label, faces_label, faces_label, elems_label]
    @test bc_label(msh_file) == [bc₁_label, bc₂_label, bc₃_label, bc₄_label, ""]

    entity_index = 100
    @test entity_label(msh_file, entity_index) == elems_label
    @test material_label(msh_file, entity_index) == "svkHyper"
    @test bc_label(msh_file, entity_index) == ""
    @test physical_index(msh_file, entity_index) == 5

    # Create the mesh
    tetra_type = [Tetrahedron(elems_label)]
    triangle_face_type = [TriangularFace(faces_label)]
    entities = StructuralEntity(tetra_type, triangle_face_type)
    gmsh_mesh = Mesh(msh_file, entities)

    # Replace first node
    new_coordinates = Point(1.0, 1.0, 1.0)
    node_index_to_replace = 1
    replace!(gmsh_mesh, node_index_to_replace, new_coordinates)
    @test coordinates(node(gmsh_mesh, node_index_to_replace)) == new_coordinates
end

@testset "ONSAS.Meshes Sets" begin
    Lᵢ = rand() * 20
    Lⱼ = rand() * 20
    Lₖ = rand() * 20

    # Mesh
    n₁ = Node(0.0, 0.0, 0.0)
    n₂ = Node(0.0, 0.0, Lₖ)
    n₃ = Node(0.0, Lⱼ, Lₖ)
    n₄ = Node(0.0, Lⱼ, 0.0)
    n₅ = Node(Lᵢ, 0.0, 0.0)
    n₆ = Node(Lᵢ, 0.0, Lₖ)
    n₇ = Node(Lᵢ, Lⱼ, Lₖ)
    n₈ = Node(Lᵢ, Lⱼ, 0.0)
    vec_nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
    # nothing is a placeholder for extra data
    mesh = Mesh(; nodes = vec_nodes)
    ## Faces
    f₁ = TriangularFace(n₅, n₈, n₆)
    f₂ = TriangularFace(n₆, n₈, n₇)
    f₃ = TriangularFace(n₄, n₁, n₂)
    f₄ = TriangularFace(n₄, n₂, n₃)
    f₅ = TriangularFace(n₆, n₂, n₁)
    f₆ = TriangularFace(n₆, n₁, n₅)
    f₇ = TriangularFace(n₁, n₄, n₅)
    f₈ = TriangularFace(n₄, n₈, n₅)
    vec_faces = [f₁, f₂, f₃, f₄, f₅, f₆, f₇, f₈]
    append!(faces(mesh), vec_faces)
    ## Entities
    t₁ = Tetrahedron(n₁, n₄, n₂, n₆)
    t₂ = Tetrahedron(n₆, n₂, n₃, n₄)
    t₃ = Tetrahedron(n₄, n₃, n₆, n₇)
    t₄ = Tetrahedron(n₄, n₁, n₅, n₆)
    t₅ = Tetrahedron(n₄, n₆, n₅, n₈)
    t₆ = Tetrahedron(n₄, n₇, n₆, n₈)
    vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
    append!(elements(mesh), vec_elems)

    # Sets
    # nodes
    node_ids_in_vec = [1, 2, 3, 4]
    nodes_set = nodes(mesh)[node_ids_in_vec]
    node_set_label = "left"
    # add using node indexes
    [add_node_to_set!(mesh, node_set_label, i) for i in node_ids_in_vec[1:3]]
    # add using the node itself
    add_node_to_set!(mesh, node_set_label, last(nodes_set))
    @test all([i ∈ node_set(mesh, node_set_label) for i in node_ids_in_vec])
    @test all([n ∈ nodes(mesh, node_set_label) for n in nodes_set])

    # faces
    face_ids_in_vec = [6, 5]
    face_set_label = "front"
    faces_set = faces(mesh)[face_ids_in_vec]
    # add using face indexes
    add_face_to_set!(mesh, face_set_label, first(face_ids_in_vec))
    # add using the face itself
    add_face_to_set!(mesh, face_set_label, last(faces_set))
    @test all([i ∈ face_set(mesh, face_set_label) for i in face_ids_in_vec])
    @test all([n ∈ faces(mesh, face_set_label) for n in faces_set])

    elem_ids_in_vec = [5, 6]
    element_set_label = "elems-with-node8"
    elements_set = elements(mesh)[elem_ids_in_vec]
    # add using element indexes
    add_element_to_set!(mesh, element_set_label, first(elem_ids_in_vec))
    # add using the element itself
    add_element_to_set!(mesh, element_set_label, last(elements_set))
    @test all([i ∈ element_set(mesh, element_set_label) for i in elem_ids_in_vec])
    @test all([n ∈ elements(mesh, element_set_label) for n in elements_set])
end
