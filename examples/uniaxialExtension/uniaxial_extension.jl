# -------------------------------------------------------------------------- 
# Uniaxial Extension ExampleExercise 4 from section 6.5 in (Holzapfel,2000).
# --------------------------------------------------------------------------
using ONSAS.StaticAnalyses
using Test: @test
## scalar parameters
# scalar parameters
E = 1.0 # Young modulus in Pa
ν = 0.3  # Poisson's ratio
p = 3   # Tension load in Pa
Lᵢ = 2.0   # Dimension in x of the box in m 
Lⱼ = 1.0   # Dimension in y of the box in m
Lₖ = 1.0   # Dimension in z of the box in m
const RTOL = 1e-4               # Relative tolerance for tests
# -------------------------------
# Mesh
#--------------------------------
n₁ = Node(0.0, 0.0, 0.0)
n₂ = Node(0.0, 0.0, Lₖ)
n₃ = Node(0.0, Lⱼ, Lₖ)
n₄ = Node(0.0, Lⱼ, 0.0)
n₅ = Node(Lᵢ, 0.0, 0.0)
n₆ = Node(Lᵢ, 0.0, Lₖ)
n₇ = Node(Lᵢ, Lⱼ, Lₖ)
n₈ = Node(Lᵢ, Lⱼ, 0.0)
vec_nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
s_mesh = Mesh(vec_nodes)
## Faces 
f₁ = TriangularFace(n₅, n₈, n₆, "loaded_face_1")
f₂ = TriangularFace(n₆, n₈, n₇, "loaded_face_2")
f₃ = TriangularFace(n₄, n₁, n₂, "x=0_face_1")
f₄ = TriangularFace(n₄, n₂, n₃, "x=0_face_2")
f₅ = TriangularFace(n₆, n₂, n₁, "y=0_face_1")
f₆ = TriangularFace(n₆, n₁, n₅, "y=0_face_2")
f₇ = TriangularFace(n₁, n₄, n₅, "z=0_face_1")
f₈ = TriangularFace(n₄, n₈, n₅, "z=0_face_2")
vec_faces = [f₁, f₂, f₃, f₄, f₅, f₆, f₇, f₈]
push!(s_mesh, vec_faces)
## Elements 
t₁ = Tetrahedron(n₁, n₄, n₂, n₆, "tetra_1")
t₂ = Tetrahedron(n₆, n₂, n₃, n₄, "tetra_2")
t₃ = Tetrahedron(n₄, n₃, n₆, n₇, "tetra_3")
t₄ = Tetrahedron(n₄, n₁, n₅, n₆, "tetra_4")
t₅ = Tetrahedron(n₄, n₆, n₅, n₈, "tetra_5")
t₆ = Tetrahedron(n₄, n₇, n₆, n₈, "tetra_6")
vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
push!(s_mesh, vec_elems)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 3
add!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
svk = SVK(E, ν, "steel")
mat_dict = dictionary([svk => [t₁, t₂, t₃, t₄, t₅, t₆]])
s_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁ = FixedDofBoundaryCondition([:u], [1], "fixed_uₓ")
bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed_uⱼ")
bc₃ = FixedDofBoundaryCondition([:u], [3], "fixed_uₖ")
# Load
bc₄ = GlobalLoadBoundaryCondition([:u], t -> [p * t, 0, 0], "tension")
# Assign this to faces 
face_bc = dictionary([bc₁ => [f₃, f₄], bc₂ => [f₅, f₆], bc₃ => [f₇, f₈], bc₄ => [f₁, f₂]])
# Crete boundary conditions struct
s_boundary_conditions = StructuralBoundaryConditions(face_bcs=face_bc)
# -------------------------------
# Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_boundary_conditions)
# -------------------------------
# Structural Analysis
# -------------------------------
# Final load factor
NSTEPS = 8
sa = StaticAnalysis(s, NSTEPS=NSTEPS)
# -------------------------------
# Algorithm
# -------------------------------
tol_f = 1e-8;
tol_u = 1e-8;
max_iter = 30;
tols = ConvergenceSettings(tol_u, tol_f, max_iter)
nr = NewtonRaphson(tols)
# -------------------------------
# Numerical solution
# -------------------------------
states_sol = solve(sa, nr)
numerical_uₖ = getindex.(displacements.(states_sol), index(Dof(19)))
numerical_λᵥ = load_factors(sa)
#-----------------------------
# Analytic solution  
#-----------------------------
"Analytic load factor solution for the displacement `uᵢ` towards `x` axis at node `n₆`."
load_factors_analytic(uᵢ::Real, p::Real=p, E::Real=E, Lᵢ::Real=Lᵢ) =
    1 / p * E * 0.5 * ((1 + uᵢ / Lᵢ)^3 - (1 + uᵢ / Lᵢ))
analytics_λᵥ = load_factors_analytic.(numerical_uₖ)
#-----------------------------
# Test boolean for CI  
#-----------------------------
@test analytics_λᵥ ≈ numerical_λᵥ rtol = RTOL