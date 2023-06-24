using Test
using Dictionaries

using ONSAS.Handlers
using ONSAS.Interpolators
using ONSAS.Entities
using ONSAS.Nodes
using ONSAS.TriangularFaces
using ONSAS.Tetrahedrons
using ONSAS.Meshes
using ONSAS.Searches

const RTOL = 1e-5

@testset "ONSAS.Meshes.PointEvalHandler" begin
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
    mesh = Mesh(; nodes=vec_nodes)
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

    # Dofs
    dof_dim = 3
    apply!(mesh, :u, dof_dim)

    # Interpolator
    #--------------------------------
    # Create a node otside the mesh
    n₉ = Node((n₇ + n₈)...)
    # Inner nodes
    nodes_to_interpolate = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈, n₉]
    vec_points = coordinates.(nodes_to_interpolate)

    # Create point eval handler and tests
    ph_nodes = PointEvalHandler(mesh, vec_points; alg=Partition())

    @test ph_nodes.mesh == mesh
    in_mesh_indexes = [1, 2, 3, 4, 5, 6, 7, 8]
    @test points(ph_nodes) == view(nodes_to_interpolate, in_mesh_indexes)
    @test n₉ ∈ not_in_mesh_points(ph_nodes)
    point_elements = points_to_element(interpolator(ph_nodes))
    node_2_weights = node_to_weights(interpolator(ph_nodes))
    @test length(point_elements) == length(points(ph_nodes)) == length(in_mesh_indexes)

    # If the point is at two elements then the first element will be reported
    @test all([p ∈ point_elements[index_p] for (index_p, p) in enumerate(points(ph_nodes))])

    # Test interpolation for a linear scalar field
    linear_scalar_field(x, y, z) = 10x - 2y + 3z + 12
    vec_linear_scalar_field = dictionary([n => linear_scalar_field(coordinates(n)...)
                                          for n in nodes_to_interpolate])

    p₉ = Point(rand() * Lᵢ, rand() * Lⱼ, rand() * Lₖ)
    ph_rand = PointEvalHandler(mesh, p₉)
    interpol = interpolator(ph_rand)
    node_2_weights = node_to_weights(interpol)
    node_to_interpolate_p₉ = first(node_2_weights)
    node_value_weighted = [weight * vec_linear_scalar_field[node]
                           for (node, weight) in pairs(node_to_interpolate_p₉)]
    @test linear_scalar_field(p₉...) ≈ reduce(+, node_value_weighted) rtol = RTOL
end
