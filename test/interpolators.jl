using Test, Dictionaries
using ONSAS.Interpolators
using ONSAS.Entities
using ONSAS.Nodes
using ONSAS.Meshes

@testset "ONSAS.Interpolators.FEMInterpolator" begin
    Lᵢ = 2.0  # Dimension in x of the box in m 
    Lⱼ = 1.0  # Dimension in y of the box in m
    Lₖ = 1.0  # Dimension in z of the box in m

    # Mesh
    ## Nodes
    n₁ = Node(0.0, 0.0, 0.0)
    n₂ = Node(0.0, 0.0, Lₖ)
    n₃ = Node(0.0, Lⱼ, Lₖ)
    n₄ = Node(0.0, Lⱼ, 0.0)
    n₅ = Node(Lᵢ, 0.0, 0.0)
    n₆ = Node(Lᵢ, 0.0, Lₖ)
    n₇ = Node(Lᵢ, Lⱼ, Lₖ)
    n₈ = Node(Lᵢ, Lⱼ, 0.0)
    nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
    ## Faces 
    f₁ = TriangularFace(n₅, n₈, n₆)
    f₂ = TriangularFace(n₆, n₈, n₇)
    f₃ = TriangularFace(n₄, n₁, n₂)
    f₄ = TriangularFace(n₄, n₂, n₃)
    f₅ = TriangularFace(n₆, n₂, n₁)
    f₆ = TriangularFace(n₆, n₁, n₅)
    f₇ = TriangularFace(n₁, n₄, n₅)
    f₈ = TriangularFace(n₄, n₈, n₅)
    faces = [f₁, f₂, f₃, f₄, f₅, f₆, f₇, f₈]
    ## Entities 
    t₁ = Tetrahedron(n₁, n₄, n₂, n₆)
    t₂ = Tetrahedron(n₆, n₂, n₃, n₄)
    t₃ = Tetrahedron(n₄, n₃, n₆, n₇)
    t₄ = Tetrahedron(n₄, n₁, n₅, n₆)
    t₅ = Tetrahedron(n₄, n₆, n₅, n₈)
    t₆ = Tetrahedron(n₄, n₇, n₆, n₈)
    elements = [t₁, t₂, t₃, t₄, t₅, t₆]

    ## Mesh
    mesh = Mesh(; nodes, elements, faces)

    ## points to interpolate 
    p₁ = Point(coordinates(n₁)...)
    p₂ = Point(coordinates(n₃)...)
    vec_points = [p₁, p₂]

    ## weights to interpolate each point
    node_to_w = [dictionary([n₁ => 1.0, n₄ => 0.0, n₂ => 0.0, n₆ => 0.0]),
                 dictionary([n₆ => 0.0, n₂ => 0.0, n₃ => 1.0, n₄ => 0.0])]

    ## Entities where each point is located
    point_to_elem = [t₁, t₂]

    ## Constructor tests
    fem_interpolator = FEMInterpolator(vec_points, node_to_w, point_to_elem)
    @test points(fem_interpolator) == vec_points
    @test node_to_weights(fem_interpolator) == node_to_w
    @test points_to_element(fem_interpolator) == point_to_elem

    ## Interpolate a scalar magnitude
    nodal_scalar_magnitude = dictionary([n₁ => 1.0, n₂ => 2.0, n₃ => 3.0,
                                         n₄ => 4.0, n₅ => 5.0, n₆ => 6.0,
                                         n₇ => 7.0, n₈ => 8.0])

    ## Manufactured interpolation
    interpolated_magnitude = zeros(length(vec_points))
    for point_idx in 1:length(points(fem_interpolator))
        for (node, weight) in pairs(node_to_w[point_idx])
            interpolated_magnitude[point_idx] += weight * nodal_scalar_magnitude[node]
        end
    end
    @test interpolate(nodal_scalar_magnitude, fem_interpolator) == interpolated_magnitude

    ## Interpolate a dim-dimensional magnitude
    nodal_magnitude = dictionary([n₁ => [1.0, 0.0, 0.0],
                                  n₂ => [2.0, 0.0, 0.0],
                                  n₃ => [3.0, 0.0, 0.0],
                                  n₄ => [4.0, 0.0, 0.0],
                                  n₅ => [5.0, 0.0, 0.0],
                                  n₆ => [6.0, 0.0, 0.0],
                                  n₇ => [7.0, 0.0, 0.0],
                                  n₈ => [8.0, 0.0, 0.0]])

    # Manufactured interpolation
    nodal_mag_dim = length(first(values(nodal_magnitude)))
    num_points_to_interpolate = length(point_to_elem)
    interpolated_magnitude = [zeros(nodal_mag_dim) for _ in 1:num_points_to_interpolate]

    for point_idx in 1:length(points(fem_interpolator))
        for (node, weight) in pairs(node_to_w[point_idx])
            interpolated_magnitude[point_idx] .+= weight * nodal_magnitude[node]
        end
    end

    @test interpolate(nodal_magnitude, fem_interpolator) == interpolated_magnitude
end
