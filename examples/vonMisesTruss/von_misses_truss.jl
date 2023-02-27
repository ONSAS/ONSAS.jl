# ------------------------------------------------------------- 
# Von Misses Truss Example from (Zerpa, Bazzano 2017 ) - 2.5.4
# -------------------------------------------------------------
using Test: @test
using ONSAS.StaticAnalyses
## scalar parameters
const E = 210e9                 # Young modulus in Pa
const ν = 0.0                   # Poisson's modulus
const A₀ = 2.5e-3                # Cross-section area in m²
const ANG = 65                  # truss angle in degrees
const L = 2                     # Length in m 
const V = L * cos(deg2rad(ANG)) # vertical distance in m 
const H = L * sin(deg2rad(ANG)) # horizontal distance in m
const Fₖ = -3e8                 # Vertical load in N
const RTOL = 1e-4               # Relative tolerance for tests
# -------------
# Mesh
# -------------
## Nodes
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(V, 0.0, H)
n₃ = Node(2V, 0.0, 0.0)
vec_nodes = [n₁, n₂, n₃]
## Cross sections
d = sqrt(4 * A₀ / pi)
s₁ = Circle(d)
a = sqrt(A₀)
s₂ = Square(a)
## Elements 
truss₁ = Truss(n₁, n₂, s₁, "left_truss") # [n₁, n₂]
truss₂ = Truss(n₂, n₃, s₂, "right_truss") # [n₂, n₃]
vec_elems = [truss₁, truss₂]
## Mesh
s_mesh = Mesh(vec_nodes, vec_elems)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 3
add!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν, "steel")
mat_dict = dictionary([steel => [truss₁, truss₂]])
s_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁ = FixedDofBoundaryCondition([:u], [1, 2, 3], "fixed_uₓ_uⱼ_uₖ")
bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed_uⱼ")
# Load 
bc₃ = GlobalLoadBoundaryCondition([:u], t -> [0, 0, Fₖ * t], "load in j")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂]])
s_boundary_conditions = StructuralBoundaryConditions(node_bcs=node_bc)
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
# Algorithm
# -------------------------------
tol_f = 1e-7;
tol_u = 1e-7;
max_iter = 100;
tols = ConvergenceSettings(tol_u, tol_f, max_iter)
nr = NewtonRaphson(tols)
# -------------------------------
# Numerical solution
# -------------------------------
states_sol = solve(sa, nr)
numerical_uₖ = getindex.(displacements.(states_sol), index(Dof(6)))
numerical_λᵥ = -load_factors(sa) * Fₖ
#-----------------------------
# Analytic solution  
#-----------------------------
"Analytic load factor solution for the displacement `uₖ` towards z axis at node `n₂`."
function load_factors_analytic(uₖ::Real, E::Real=E, A::Real=A₀, H::Real=H, V::Real=V, l₀=L)
    λ = -2 * E * A *
        ((H + uₖ)^2 + V^2 - l₀^2) /
        (l₀ * (l₀ + sqrt((H + uₖ)^2 + V^2))) *
        (H + uₖ) / sqrt((H + uₖ)^2 + V^2)
end
analytics_λᵥ = load_factors_analytic.(numerical_uₖ)
#-----------------------------
# Test boolean for CI  
#-----------------------------
@test analytics_λᵥ ≈ numerical_λᵥ rtol = RTOL
