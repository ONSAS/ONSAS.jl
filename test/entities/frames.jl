using Test
using ONSAS.Nodes
using ONSAS.Rectangles
using ONSAS.Frames

@testset "Frame elements." begin
    L1 = 2.0
    L2 = 1.5
    ty = 0.1
    tz = 0.2

    n1 = Node(0.0, L1, L1)
    n2 = Node(0.0, 0.0, L2)
    n3 = Node(0.0, 0.0, 0.0)
    R = Rectangle(ty, tz)
    f1 = Frame(n1, n2, R)

    @test nodes(f1) == [n1, n2]
    @test cross_section(f1) == R
    @test local_dof_symbol(f1) == [:u, :θ]
    @test f1.mass_matrix == Consistent

    f2 = Frame(n1, n2, R, Lumped)
    @test f2.mass_matrix == Lumped

    # Constructor with different types
    n1 = Node(1.2, 2.3, 3.4)
    n2 = Node(2, 3, 4)
    S = Rectangle(0.6, 0.3)
    f = Frame(n1, n2, S)
    @test coordinates.(nodes(f)) == coordinates.([n1, n2])
end
