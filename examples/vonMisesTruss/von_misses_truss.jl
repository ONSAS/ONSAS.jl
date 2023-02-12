## Von Mises truss example problem
using ONSAS
using LinearAlgebra
using Test
## scalar parameters
E = 2e11  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 5e-3  # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m 
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
# Fx = 0     # horizontal load in N
Fⱼ = -3e8  # vertical   load in N
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν, "steel")
aluminum = SVK(E / 3, ν, "aluminium")
# -------------------------------
# Geometries
# -------------------------------
## Cross section
a = sqrt(A)
s₁ = Square(a)
s₂ = Square(2a)
# -------------------------------
# Create mesh
# -------------------------------
## Nodes
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(d, h, 0.0)
n₃ = Node(2d, 0.0, 0.0)
vec_nodes = [n₁, n₂, n₃]
## Elements 
truss₁ = Truss(n₁, n₂, s₁, "left_truss") # [n₁, n₂]
truss₂ = Truss(n₂, n₃, s₂, "right_truss") # [n₂, n₃]
vec_elems = [truss₁, truss₂]
## Mesh
s_mesh = Mesh(vec_nodes, vec_elems)
# -------------------------------
# Dofs
#--------------------------------}
# add!(s_mesh, :u, 3)




#=

# -------------------------------
# Materials
# -------------------------------
s_materials = StructuralMaterials(
    dictionary([steel => [truss₁], aluminum => [truss₂]])
)
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = PinnedDisplacementBoundaryCondition("fixed")
bc₂ = FⱼLoadBoundaryCondition(Fⱼ, "load in j")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂]])
elem_bc = dictionary([bc₁ => [truss₁], bc₂ => [truss₂]])
s_boundary_conditions = StructuralBoundaryConditions(node_bc, elem_bc)

for (bc, n) in pairs(node_bc)
    apply!(s, bc, n)
end
# -------------------------------
# Create Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_boundary_conditions)

# -------------------------------
# Structural Analysis
# -------------------------------
# Final load factor
λ₁ = 10
NSTEPS = 9
s_analysis = StaticAnalysis(s, λ₁, NSTEPS=NSTEPS)
s_state = current_state(s_analysis)
@test norm(displacements(s_state)) == 0
@test norm(external_forces(s_state)) == 0
@test norm(internal_forces(s_state)) == 0
# -------------------------------
# Solve analysis
# -------------------------------
tol_f = 1e-10;
tol_u = 1e-10;
max_iter = 100;
tols = ConvergenceSettings(tol_f, tol_u, max_iter)
nr = NewtonRaphson(tols)
sol = solve(s_analysis, nr)
# typeof(sol) = StaticSolution(:u, hist:)

=#