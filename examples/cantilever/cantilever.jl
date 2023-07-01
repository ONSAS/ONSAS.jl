using Test
using ONSAS

# ------------------------
# Model parameters
# ------------------------
b = 0.2 # Width
h = 0.3 # Height
pp = 2400
q = b * h * pp * 9.81 # Distributed load
Px = 1e2
Py = 1e3
E = 210e9 # Elastic modulus
RTOL = 1e-4 # Relative tolerance for tests

# ------------------------
# Mesh
# ------------------------
L = 3
num_elems = 5
x_coords = range(0, L, num_elems + 1)
nodes = [Node(xi, 0.0, 0.0) for xi in x_coords]
S = Rectangle(h, b)
frames = [Frame(nodes[j], nodes[j + 1], S) for j in 1:(length(nodes) - 1)]
msh = Mesh(; nodes=nodes, elements=frames)
set_dofs!(msh, :u, 3)
set_dofs!(msh, :θ, 3)

# ------------------------
# Materials and BC
# ------------------------
i = IsotropicLinearElastic(E, 0.3)
mat = StructuralMaterial(i => frames)
bc1 = FixedDof(:u, [1, 2, 3])
bc2 = FixedDof(:θ, [1, 2, 3])
bc3 = GlobalLoad(:u, t -> [Px, -Py, 0])
bc = StructuralBoundaryCondition(bc1 => [nodes[1]], bc2 => [nodes[1]], bc3 => [nodes[end]])

# ------------------------
# Perform linear analysis
# ------------------------
s = Structure(msh, mat, bc)
anali = LinearStaticAnalysis(s; NSTEPS=10)
sol = solve!(anali)

# ------------------------
# Verifications
# ------------------------

@testset "Displacements at the right-most node." begin
    Izz = b * h^3 / 12
    A = b * h

    @test displacements(sol, nodes[end], 2)[1] ≈ -Py * L^3 / (3 * E * Izz) rtol = RTOL
    @test displacements(sol, nodes[end], 6)[1] ≈ -Py * L^2 / (2 * E * Izz) rtol = RTOL
    @test displacements(sol, nodes[end], 1)[1] ≈ Px * L / (E * A) rtol = RTOL
end
