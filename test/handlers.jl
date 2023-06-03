using Test
using ONSAS.Handlers
using ONSAS.Elements
using ONSAS.Meshes
using ONSAS.Interpolators

const RTOL = 1e-5

@testset "ONSAS.Meshes.PointEvalHandler + TriangularFace + Tetrahedron + Sets" begin
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
    ## Elements 
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

    # Dofs
    dof_dim = 3
    apply!(mesh, :u, dof_dim)

    # Interpolator
    #--------------------------------
    # Outter nodes
    n₉ = Node((n₇ + n₈)...)
    # Inner nodes
    nodes_to_interpolate = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈, n₉]
    vec_points = coordinates.(nodes_to_interpolate)
    ph_nodes = PointEvalHandler(mesh, vec_points)

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
