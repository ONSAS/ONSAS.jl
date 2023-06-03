#######################
# Meshes module tests #
#######################
using Test, Suppressor, Dictionaries
using ONSAS.Meshes
using ONSAS.Gmsh
using ONSAS.Elements
using ONSAS.Circles

const RTOL = 1e-5

@testset "ONSAS.Meshes.Mesh Node + Trusses" begin

    # Mesh considering only nodes
    # Nodes
    L = 1
    n₁ = Node(0, 0, 0)
    n₂ = Node(L, L, L)
    n₃ = Node(2L, 0, 5L)
    n₄ = Node(2L, 0, 5L)
    vec_nodes = [n₁, n₂, n₃, n₄]

    ## Elements
    d = 0.2
    s₁ = Circle(d)
    t₁ = Truss(n₁, n₂, s₁)
    t₂ = Truss(n₂, n₃, s₁)
    t₃ = Truss(n₃, n₄, s₁)
    vec_elements = [t₁, t₂, t₃]

    # Constructors
    mesh = Mesh(; nodes=vec_nodes, elements=vec_elements)

    @test dimension(mesh) == dimension(n₁)
    @test all(isempty.(dofs(mesh)))
    @test elements(mesh) == vec_elements
    @test nodes(mesh) == vec_nodes

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
    apply!(mesh, :u, u_dim)
    @test num_nodes(mesh) * u_dim == num_dofs(mesh)
    θ_dim = 3
    apply!(mesh, :θ, θ_dim)
    @test num_nodes(mesh) * (u_dim + θ_dim) == num_dofs(mesh)
end

uniaxial_mesh_path = joinpath(@__DIR__, "..", "examples", "uniaxial_extension", "uniaxial_mesh.jl")
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
end
