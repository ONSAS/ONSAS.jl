using Test
using StaticArrays

using ONSAS.Searches
using ONSAS.Entities
using ONSAS.Nodes
using ONSAS.Meshes
using ONSAS.TriangularFaces
using ONSAS.Tetrahedrons

@testset "ONSAS.Searches" begin
    # -------------------------------
    # Mesh
    #--------------------------------
    Lᵢ = 2.0
    Lⱼ = 1.0
    Lₖ = 1.0
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
    f₁ = TriangularFace(n₅, n₈, n₆, "loaded_face_1")
    f₂ = TriangularFace(n₆, n₈, n₇, "loaded_face_2")
    f₃ = TriangularFace(n₄, n₁, n₂, "x=0_face_1")
    f₄ = TriangularFace(n₄, n₂, n₃, "x=0_face_2")
    f₅ = TriangularFace(n₆, n₂, n₁, "y=0_face_1")
    f₆ = TriangularFace(n₆, n₁, n₅, "y=0_face_2")
    f₇ = TriangularFace(n₁, n₄, n₅, "z=0_face_1")
    f₈ = TriangularFace(n₄, n₈, n₅, "z=0_face_2")
    faces = [f₁, f₂, f₃, f₄, f₅, f₆, f₇, f₈]
    ## Elements
    t₁ = Tetrahedron(n₁, n₄, n₂, n₆, "tetra_1")
    t₂ = Tetrahedron(n₆, n₂, n₃, n₄, "tetra_2")
    t₃ = Tetrahedron(n₄, n₃, n₆, n₇, "tetra_3")
    t₄ = Tetrahedron(n₄, n₁, n₅, n₆, "tetra_4")
    t₅ = Tetrahedron(n₄, n₆, n₅, n₈, "tetra_5")
    t₆ = Tetrahedron(n₄, n₇, n₆, n₈, "tetra_6")
    elements = [t₁, t₂, t₃, t₄, t₅, t₆]
    mesh = Mesh(; nodes, elements, faces)

    # -------------------------------
    # Search
    #--------------------------------
    points_to_check = [Point(coordinates(n₁)...),
        Point(n₁ .- SVector(1, 1, 1)),
        Point(n₇ .+ SVector(10, 10, 10)),
        Point(coordinates(n₇)...),
        Point(coordinates(n₄)...)]

    # Serial implementation
    in_mesh_points_idx,
    in_mesh_elements_idx = evaluate_points_in_mesh(mesh,
        points_to_check,
        Serial())

    @test all([p in mesh for p in view(points_to_check, in_mesh_points_idx)])
    @test in_mesh_elements_idx == [1, 3, 1]

    # Parallel implementation
    in_mesh_points_idx,
    in_mesh_elements_idx = evaluate_points_in_mesh(mesh,
        points_to_check,
        Threaded())
    @test all([p in mesh for p in view(points_to_check, in_mesh_points_idx)])
    # the order can be altered
    @test all(index in in_mesh_elements_idx for index in [1, 3, 1])

    # Partitioned implementation
    in_mesh_points_idx,
    in_mesh_elements_idx = evaluate_points_in_mesh(mesh,
        points_to_check,
        Partition())
    @test all([p in mesh for p in view(points_to_check, in_mesh_points_idx)])
    @test in_mesh_elements_idx == [1, 3, 1]

    # Partitioned thread implementation
    in_mesh_points_idx,
    in_mesh_elements_idx = evaluate_points_in_mesh(mesh,
        points_to_check,
        PartitionThreaded())
    @test all([p in mesh for p in view(points_to_check, in_mesh_points_idx)])
    # the order can be altered
    @test all(index in in_mesh_elements_idx for index in [1, 3, 1])
end
