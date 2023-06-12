using Test
using Suppressor
using ONSAS.Entities
using ONSAS.Nodes
using ONSAS.TriangularFaces
using ONSAS.Tetrahedrons
using ONSAS.Gmsh
using ONSAS.Meshes
using ONSAS.StructuralEntities

path = joinpath("..", "..", "examples", "uniaxial_extension", "uniaxial_mesh.jl")
include(path)

@testset "ONSAS.Meshes.GMSH.MshFile " begin

    # MshFile
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
    rm(file_name; force=true)

    @test nodes(msh_file) == msh_file.vec_nodes
    @test length(physical_index(msh_file)) == length(connectivity(msh_file))
    @test dimension(msh_file) == dimension(first(nodes(msh_file)))
    @test material_label(msh_file) == ["", "", "", "", "svkHyper"]
    @test entity_label(msh_file) == ["triangle", "triangle", "triangle", "triangle", "tetrahedron"]
    @test bc_label(msh_file) == ["fixed-ux", "fixed-uj", "fixed-uk", "tension", ""]

    entity_index = 100
    @test entity_label(msh_file, entity_index) == "tetrahedron"
    @test material_label(msh_file, entity_index) == "svkHyper"
    @test bc_label(msh_file, entity_index) == ""
    @test physical_index(msh_file, entity_index) == 5

    # Mesh
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elems_label)]
    structural_entities = StructuralEntity(velems, vfaces)

    msh_mesh = Mesh(msh_file, structural_entities)

    @test Meshes.nodes(msh_mesh) === nodes(msh_file)

    # pick an element corresponding to the  samae type as entity_index
    element_index = element_set(msh_mesh, String(entity_label(msh_file, entity_index)))
    test_element = element(msh_mesh, rand(element_index))
    @test String(label(test_element)) == entity_label(msh_file, entity_index)

    # Replace node
    node_index_in_the_element = 1
    old_node_element = nodes(test_element)[node_index_in_the_element]
    old_node_element ∈ nodes(msh_mesh)

    node_index = findfirst(x -> x === old_node_element, nodes(msh_mesh))
    new_coordinates = Point(1.0, 1.0, 1.0)
    new_node = Node(new_coordinates, dofs(old_node_element))
    Meshes.nodes(msh_mesh)[node_index] = new_node
    @test nodes(msh_mesh)[node_index] === new_node

    @test nodes(test_element)[node_index_in_the_element] === new_node
end
