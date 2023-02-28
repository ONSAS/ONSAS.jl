# -------------------------------------------------------------------------- 
# Uniaxial Extension ExampleExercise 4 from section 6.5 in (Holzapfel,2000).
# For notation see: https://onsas.github.io/ONSAS.m/dev/examples/uniaxialExtension/
# --------------------------------------------------------------------------
using ONSAS.StaticAnalyses
using ONSAS.Utils: eye
using Test: @test
using LinearAlgebra: det, tr
using Roots: find_zero
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
# Displacements in the x (component 1) axis at node 7
numerical_uᵢ = getindex.(displacements.(states_sol), index(dofs(n₇)[:u][1]))
numerical_α = 1 .+ numerical_uᵢ / Lᵢ
# Displacements in the y (component 2) axis at node 7
numerical_uⱼ = getindex.(displacements.(states_sol), index(dofs(n₇)[:u][2]))
numerical_β = 1 .+ numerical_uⱼ / Lⱼ
# Displacements in the z (component 3) axis at node 7
numerical_uₖ = getindex.(displacements.(states_sol), index(dofs(n₇)[:u][3]))
numerical_γ = 1 .+ numerical_uᵢ / Lₖ
# Extract ℙ and ℂ from the last state
element_index = 5
# Cosserat or second Piola-Kirchhoff stress tensor
ℙ_numeric = collect(values(stress(last(states_sol))))[element_index]
# Right hand Cauchy strain tensor 
ℂ_numeric = collect(values(strain(last(states_sol))))[element_index]
# Load factors 
numerical_λᵥ = load_factors(sa)
#-----------------------------
# Analytic solution  
#-----------------------------
# Test with load factors
"Analytic load factor solution for the displacement `uᵢ` towards `x` axis at node `n₆`."
load_factors_analytic(uᵢ::Real, p::Real=p, E::Real=E, Lᵢ::Real=Lᵢ) =
    1 / p * E * 0.5 * ((1 + uᵢ / Lᵢ)^3 - (1 + uᵢ / Lᵢ))
analytics_λᵥ = load_factors_analytic.(numerical_uᵢ)
# Test last step σ and ϵ
α_analytic = find_zero(α -> E / 2 * α * (α^2 - 1) - p * last(load_factors(sa)), 1e-2)
β_analytic = sqrt(-ν * (α_analytic^2 - 1) + 1)
# Gradient tensor
# 𝕦 = (αx, βy, γz)
𝔽_analytic = [
    α_analytic 0 0
    0 β_analytic 0
    0 0 β_analytic
]
# Right hand Cauchy tensor 
ℂ_analytic = 𝔽_analytic' * 𝔽_analytic
𝕁 = det(ℂ_analytic)
# Green-Lagrange strain tensor
𝕀 = eye(3)
𝔼_analytic = 1 / 2 * (ℂ_analytic - 𝕀)
# Cosserat or second Piola-Kirchhoff stress tensor
p₁, p₂ = lame_parameters(svk)
𝕊_analytic = p₁ * tr(𝔼_analytic) * eye(3) + 2 * p₂ * 𝔼_analytic
# First Piola-Kirchhoff stress tensor
ℙ_analytic = 𝔽_analytic * 𝕊_analytic
# Cauchy stress tensor
# σ = ℙ_analytic * 𝔽_analytic'
#-----------------------------
# Test boolean for CI  
#-----------------------------
@test analytics_λᵥ ≈ numerical_λᵥ rtol = RTOL
@test ℙ_analytic ≈ ℙ_analytic rtol = RTOL
@test α_analytic ≈ last(numerical_α) rtol = RTOL
@test β_analytic ≈ last(numerical_β) rtol = RTOL
@test ℂ_numeric ≈ ℂ_analytic rtol = RTOL