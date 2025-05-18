using Test, Dictionaries, StaticArrays

using ONSAS.Nodes
using ONSAS.Entities
using ONSAS.TriangularFaces

@testset "ONSAS.Entities.TriangularFace 3D" begin
    x₁ = [0, 0, 0]
    x₂ = [1, 0, 0]
    x₃ = [0, 1, 0]

    n₁ = Node(
        x₁, dictionary([:u => [Dof(1), Dof(2), Dof(3)], :θ => [Dof(13), Dof(14), Dof(15)]]))
    n₂ = Node(
        x₂, dictionary([:u => [Dof(4), Dof(5), Dof(6)], :θ => [Dof(16), Dof(17), Dof(18)]]))
    n₃ = Node(
        x₃, dictionary([:u => [Dof(7), Dof(8), Dof(9)], :θ => [Dof(19), Dof(20), Dof(21)]]))

    face_label = "my_face"
    f_empty_nodes = TriangularFace(face_label)
    @test label(f_empty_nodes) == Symbol(face_label)
    f₁ = TriangularFace(n₁, n₂, n₃, face_label)
    f₁_no_label = TriangularFace(n₁, n₂, n₃)
    f₁ = create_entity(f_empty_nodes, [n₁, n₂, n₃])

    @test all([n in nodes(f₁) for n in [n₁, n₂, n₃]])
    @test coordinates(f₁) == [coordinates(n₁), coordinates(n₂), coordinates(n₃)]
    @test dimension(f₁) == length(x₁)
    @test all([d in dofs(f₁)[:u] for d in Dof.(1:9)])
    @test all([d in dofs(f₁)[:θ] for d in Dof.(13:21)])
    @test label(f₁) == Symbol(face_label)
    @test area(f₁) == 0.5
    @test normal_direction(f₁) == [0, 0, 1]

    # create entity for gmsh
    empty_entity = TriangularFace(face_label)
    tf = create_entity(empty_entity, [n₁, n₂, n₃])
    @test all([n in nodes(tf) for n in [n₁, n₂, n₃]])
    @test coordinates(tf) == [coordinates(n₁), coordinates(n₂), coordinates(n₃)]
    @test label(empty_entity) == label(tf)
end
