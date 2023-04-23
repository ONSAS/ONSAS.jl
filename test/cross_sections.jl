###############################
# Cross-sections module tests #
###############################
using Test
using ONSAS.CrossSections

rand_dim_1 = rand()
rand_dim_2 = rand_dim_1 + rand()

const TOLERANCE = 1e-4

@testset "ONSAS.CrossSections.Square" begin
    a = rand_dim_1
    square = Square(a)
    @test area(square) == a^2
    @test Ixx(square) ≈ a / 2 * (a / 2)^3 * (16 / 3 - 3.36 * (1 - (a / 2)^4 / (12 * (a / 2)^4))) rtol = TOLERANCE
    @test Iyy(square) == Izz(square) == a^4 / 12
    @test Ixy(square) == Ixz(square) == Iyz(square) == 0
end

@testset "ONSAS.CrossSections.Circle" begin
    d = rand_dim_1
    circle = Circle(d)
    @test area(circle) == pi * d^2 / 4
    @test Ixx(circle) == pi * d^4 / 32
    @test Iyy(circle) == Izz(circle) == pi * d^4 / 64
    @test Ixy(circle) == Ixz(circle) == Iyz(circle) == 0
end

@testset "ONSAS.CrossSections.Rectangle" begin
    h = rand_dim_1
    b = rand_dim_2
    rectangle = Rectangle(h, b)
    @test area(rectangle) == h * b
    @test Iyy(rectangle) == h * b^3 / 12
    @test Izz(rectangle) == h^3 * b / 12
    @test Ixy(rectangle) == Ixz(rectangle) == Iyz(rectangle) == 0
    a = 1 / 2 * maximum([h, b])
    b = 1 / 2 * minimum([h, b])
    @test Ixx(rectangle) == a * b^3 * (16 / 3 - 3.36 * b / a * (1 - b^4 / (12 * a^4)))
end
