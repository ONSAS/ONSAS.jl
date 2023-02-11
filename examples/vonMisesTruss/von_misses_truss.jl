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
# Create mesh
# -------------------------------
## Nodes
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(d, h, 0.0)
n₃ = Node(2d, 0.0, 0.0)
vec_nodes = [n₁, n₂, n₃]
## Elements connectivity
elem₁_nodes = [n₁, n₂]
elem₂_nodes = [n₂, n₃]
Truss(n₁, n₂)
Tetra(n₁, n₂)
vec_conec_elems = [elem₁_nodes, elem₂_nodes]
s_mesh = Mesh(vec_nodes, vec_conec_elems)
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν)
aluminum = SVK(E / 3, ν)
materials_dict = Dict("steel" => steel, "aluminum" => aluminum)
s_materials = StructuralMaterials(materials_dict)
# Sets
m_sets = Dict(
    "steel" => Set([ElementIndex(1)]),
    "aluminum" => Set([ElementIndex(2)]),
)
add_set!(s_materials, m_sets)
# -------------------------------
# Geometries
# -------------------------------
## Cross section
a = sqrt(A)
s₁ = Square(a)
s₂ = Square(2a)
# -------------------------------
# Elements
# -------------------------------
elements_dict = Dict{String,AbstractElement}(
    "truss_s₁" => Truss(s₁),
    "truss_s₂" => Truss(s₂)
)
s_elements = StructuralElements(elements_dict)
# Sets
e_sets = Dict(
    "truss_s₁" => Set([ElementIndex(1)]),
    "truss_s₂" => Set([ElementIndex(2)]),
)
add_set!(s_elements, e_sets)
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = PinnedDisplacementBoundaryCondition()
bc₂ = FⱼLoadBoundaryCondition(Fⱼ)
bcs_dict = Dict(
    "fixed" => bc₁,
    "load" => bc₂
)
s_boundary_conditions = StructuralBoundaryConditions(bcs_dict)
bc_sets = Dict(
    "fixed" => Set([NodeIndex(1), NodeIndex(3)]),
    "load" => Set([NodeIndex(2)])
)
add_set!(s_boundary_conditions, bc_sets)
# -------------------------------
# Create Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_elements, s_boundary_conditions)
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
#
# extras::Dict{Any,Any}()
#
# plot(sol)
# stress_field(sol)
#
# typeof(sol) = StaticSolution(:u, hist:)

