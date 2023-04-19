#######################
# Meshes module tests #
#######################
using Test: @testset, @test
using ONSAS.Meshes
using Dictionaries: dictionary
using StaticArrays: SVector

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
    mesh_with_sets = Mesh(vec_nodes, vec_elements)
    mesh = Mesh(vec_nodes, vec_elements)

    @test dimension(mesh) == dimension(n₁)
    @test all(isempty.(dofs(mesh)))
    @test elements(mesh) == vec_elements
    @test nodes(mesh) == vec_nodes

    # Add nodes and elements
    new_node₁ = Node(3L, 0, 5L)
    new_node₂ = Node(3L, 0, 5L)
    push!(mesh, [new_node₁, new_node₂])
    @test last(nodes(mesh)) == new_node₂
    new_t₄ = Truss(n₃, n₄, s₁)
    new_t₅ = Truss(n₂, n₁, s₁)
    push!(mesh, [new_t₄, new_t₅])
    @test last(elements(mesh)) == new_t₅
    # Add dofs 
    u_dim = 3
    add!(mesh, :u, u_dim)
    @test num_nodes(mesh) * u_dim == num_dofs(mesh)
    θ_dim = 3
    add!(mesh, :θ, θ_dim)
    @test num_nodes(mesh) * (u_dim + θ_dim) == num_dofs(mesh)

end

uniaxial_mesh_path = joinpath(@__DIR__, "..", "examples", "uniaxial_extension", "uniaxial_mesh.jl")
include(uniaxial_mesh_path)
@testset "ONSAS.Meshes.GMSH " begin


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
    filename = "uniaxial_test_mesh"
    # Create, load  and delete the mesh
    dir = @__DIR__
    # Refinement mesh factor
    ms = 0.5
    msh_file = create_uniaxial_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename, ms, dir)
    file_name = joinpath(dir, msh_file)
    msh_file = MshFile(file_name)
    rm(file_name, force=true)

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

@testset "ONSAS.Meshes.PointEvalHandler + TriangularFace + Tetrahedron" begin

    Lᵢ = rand() * 20
    Lⱼ = rand() * 20
    Lₖ = rand() * 20

    # -------------------------------
    # Mesh
    #--------------------------------
    n₁ = Node(0.0, 0.0, 0.0)
    n₂ = Node(0.0, 0.0, Lₖ)
    n₃ = Node(0.0, Lⱼ, Lₖ)
    n₄ = Node(0.0, Lⱼ, 0.0)
    n₅ = Node(Lᵢ, 0.0, 0.0)
    n₆ = Node(Lᵢ, 0.0, Lₖ)
    n₇ = Node(Lᵢ, Lⱼ, Lₖ)
    n₈ = Node(Lᵢ, Lⱼ, 0.0)
    vec_nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
    s_mesh = Mesh(vec_nodes)
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
    push!(s_mesh, vec_faces)
    ## Elements 
    t₁ = Tetrahedron(n₁, n₄, n₂, n₆)
    t₂ = Tetrahedron(n₆, n₂, n₃, n₄)
    t₃ = Tetrahedron(n₄, n₃, n₆, n₇)
    t₄ = Tetrahedron(n₄, n₁, n₅, n₆)
    t₅ = Tetrahedron(n₄, n₆, n₅, n₈)
    t₆ = Tetrahedron(n₄, n₇, n₆, n₈)
    vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
    push!(s_mesh, vec_elems)
    # -------------------------------
    # Dofs
    #--------------------------------
    dof_dim = 3
    add!(s_mesh, :u, dof_dim)

    # Interpolator
    #--------------------------------
    # At the nodes
    nodes_to_interpolate = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
    vec_points = coordinates.(nodes_to_interpolate)
    ph_nodes = PointEvalHandler(s_mesh, vec_points)
    @test mesh(ph_nodes) == s_mesh
    @test points(ph_nodes) == vec_points
    point_elements = points_to_element(interpolator(ph_nodes))

    # If the point is at two elements then the first element will 
    # be reported 
    @test all([p ∈ point_elements[index_p] for (index_p, p) in enumerate(points(ph_nodes))])

    # Test interpolation for a linear scalar field 
    linear_scalar_field(x, y, z) = 10x - 2y + 3z + 12
    vec_linear_scalar_field = dictionary([n => linear_scalar_field(coordinates(n)...) for n in nodes_to_interpolate])

    # TODO Define Point(...) constructor.
    p₉ = SVector(rand() * Lᵢ, rand() * Lⱼ, rand() * Lₖ)
    ph_rand = PointEvalHandler(s_mesh, p₉)
    interpol = interpolator(ph_rand)
    node_2_weights = node_to_weights(interpol)
    node_to_interpolate_p₉ = first(node_2_weights)
    node_value_weighted = [weight * vec_linear_scalar_field[node] for (node, weight) in pairs(node_to_interpolate_p₉)]
    @test linear_scalar_field(p₉...) ≈ reduce(+, node_value_weighted) rtol = RTOL

end
