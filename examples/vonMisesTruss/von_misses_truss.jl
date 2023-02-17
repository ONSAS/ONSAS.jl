## Von Mises truss example problem
using ONSAS
using LinearAlgebra
using Test
## scalar parameters
E = 210e9  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 2.5e-3   # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m 
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
# Fx = 0     # horizontal load in N
Fₖ = -3e8  # vertical   load in N
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν, "steel")
# -------------------------------
# Geometries
# -------------------------------
## Cross section
a = sqrt(4 * A / pi)
s = Circle(a)
# -------------------------------
# Create mesh
# -------------------------------
## Nodes
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(d, 0.0, h)
n₃ = Node(2d, 0.0, 0.0)
vec_nodes = [n₁, n₂, n₃]
## Elements 
truss₁ = Truss(n₁, n₂, s, "left_truss") # [n₁, n₂]
truss₂ = Truss(n₂, n₃, s, "right_truss") # [n₂, n₃]
vec_elems = [truss₁, truss₂]
## Mesh
s_mesh = Mesh(vec_nodes, vec_elems)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 3
add_dofs!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
mat_dict = dictionary([steel => [truss₁, truss₂]])
s_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = FixedDofBoundaryCondition([:u], collect(1:dof_dim), "fixed_uₓ_uⱼ_uₖ")
bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed_uⱼ")
bc₃ = GlobalLoadBoundaryCondition([:u], t -> [0, 0, Fₖ * t], "load in j")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂]])
s_boundary_conditions = StructuralBoundaryConditions(node_bc)
# -------------------------------
# Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_boundary_conditions)
# -------------------------------
# Structural Analysis
# -------------------------------
# Final load factor
λ₁ = 1
NSTEPS = 10
sa = StaticAnalysis(s, λ₁, NSTEPS=NSTEPS)
# -------------------------------
# Solve analysis
# -------------------------------
tol_f = 1e-10;
tol_u = 1e-10;
max_iter = 100;
tols = ConvergenceSettings(tol_f, tol_u, max_iter)
nr = NewtonRaphson(tols)
sol = solve(sa, nr)
# typeof(sol) = StaticSolution(:u, hist:)
