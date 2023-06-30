using ONSAS

# parametros
b = 0.2 # width
h = 0.3 # heigth
pp = 2400
# carga distribuida
q = b * h * pp * 9.81
Px = 1e2
Py = 1e3
E = 210e9 # elastic modulus
RTOL = 1e-4                # Relative tolerance for tests

# -----------------------------------------
# importacion de malla
L = 3
num_elems = 5
x_coords = range(0, L, num_elems + 1)
# -----------------------------------------

# Mesh
nodes = [Node(xi, 0.0, 0.0) for xi in x_coords]
S = Rectangle(h, b) # cm^2
frames = [Frame(nodes[j], nodes[j + 1], S) for j in 1:(length(nodes) - 1)]
msh = Mesh(; nodes=nodes, elements=frames)
# -----------------------------------------
set_dofs!(msh, :u, 3)
set_dofs!(msh, :θ, 3)
# Materials
i = IsotropicLinearElastic(E, 0.3)
mat = StructuralMaterial(i => frames)
# Boundary conditions
bc1 = FixedDof(:u, [1, 2, 3])
bc2 = FixedDof(:θ, [1, 2, 3])
bc3 = GlobalLoad(:u, t -> [Px, -Py, 0])

bc = StructuralBoundaryCondition(bc1 => [nodes[1]], bc2 => [nodes[1]], bc3 => [nodes[end]])

# Structure
s = Structure(msh, mat, bc)
# Analysis
anali = LinearStaticAnalysis(s; NSTEPS=10)
@time sol = solve!(anali)

# =================
# verification
Izz = b * h^3 / 12
A = b * h

@show numer_sol_deltay = displacements(sol, nodes[end], 2)[1]
@show anali_sol_deltay = -Py * L^3 / (3 * E * Izz)

@show numer_sol_angle = displacements(sol, nodes[end], 6)[1]
@show anali_sol_angle = -Py * L^2 / (2 * E * Izz)

@show numer_sol_deltax = displacements(sol, nodes[end], 1)[1]
@show anali_sol_deltax = Px * L / (E * A)

using Test
@test anali_sol_deltay ≈ numer_sol_deltay rtol = RTOL
@test anali_sol_angle ≈ numer_sol_angle rtol = RTOL
@test anali_sol_deltax ≈ numer_sol_deltax rtol = RTOL
