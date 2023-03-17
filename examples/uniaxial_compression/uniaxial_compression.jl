# ----------------------------- 
# Uniaxial Compression Example
# -----------------------------
using ONSAS.StaticAnalyses
using ONSAS.Utils: eye
using Test: @test, @testset
using LinearAlgebra: Symmetric, norm, det, tr
using Roots: find_zero
## scalar parameters
E = 1.0                    # Young modulus in Pa
ν = 0.3                    # Poisson's ratio
μ = G = E / (2 * (1 + ν))  # Second Lamé parameter 
K = E / (3 * (1 - 2 * ν))  # Bulk modulus
p = 5                      # Tension load in Pa
Lᵢ = 2.0                   # Dimension in x of the box in m 
Lⱼ = 1.0                   # Dimension in y of the box in m
Lₖ = 1.0                   # Dimension in z of the box in m
const RTOL = 1e-4          # Relative tolerance for tests
const ATOL = 1e-10         # Absolute tolerance for tests
const GENERATE_MSH = false # Boolean to generate the .msh form .geo
# -----------------------------------------------------
# Case 1 - Manufactured mesh and `NeoHookean` material
#------------------------------------------------------
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
s₁_mesh = Mesh(vec_nodes)
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
push!(s₁_mesh, vec_faces)
## Elements 
t₁ = Tetrahedron(n₁, n₄, n₂, n₆, "tetra_1")
t₂ = Tetrahedron(n₆, n₂, n₃, n₄, "tetra_2")
t₃ = Tetrahedron(n₄, n₃, n₆, n₇, "tetra_3")
t₄ = Tetrahedron(n₄, n₁, n₅, n₆, "tetra_4")
t₅ = Tetrahedron(n₄, n₆, n₅, n₈, "tetra_5")
t₆ = Tetrahedron(n₄, n₇, n₆, n₈, "tetra_6")
vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
push!(s₁_mesh, vec_elems)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 3
add!(s₁_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
# Built neo hookian material with E and ν
neo_hookean = NeoHookean(K, μ, "NeoBuiltIn")
mat_dict = dictionary([neo_hookean => [t₁, t₂, t₃, t₄, t₅, t₆]])
s₁_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁ = FixedDofBoundaryCondition([:u], [1], "fixed-ux")
bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed-uj")
bc₃ = FixedDofBoundaryCondition([:u], [3], "fixed-uk")
# Load
bc₄ = GlobalLoadBoundaryCondition([:u], t -> [p * t, 0, 0], "compression")
# Assign this to faces 
face_bc = dictionary([bc₁ => [f₃, f₄], bc₂ => [f₅, f₆], bc₃ => [f₇, f₈], bc₄ => [f₁, f₂]])
# Crete boundary conditions struct
s₁_boundary_conditions = StructuralBoundaryConditions(face_bcs=face_bc)
# -------------------------------
# Structure
# -------------------------------
s₁ = Structure(s₁_mesh, s₁_materials, s₁_boundary_conditions)
# -------------------------------
# Structural Analysis
# -------------------------------
# Final load factor
NSTEPS = 8
sa₁ = StaticAnalysis(s₁, NSTEPS=NSTEPS)
# Resets the analysis in order to run it multiple times
reset!(sa₁)
# -------------------------------
# Algorithm
# -------------------------------
tol_f = 1e-10;
tol_u = 1e-10;
max_iter = 30;
tols = ConvergenceSettings(tol_u, tol_f, max_iter)
nr = NewtonRaphson(tols)
# -------------------------------
# Numerical solution
# -------------------------------
states_sol_case₁ = solve(sa₁, nr)
"Computes numeric solution α(L_def/L_ref), β(L_def/L_ref) and γ(L_def/L_ref)
 for analytic validation."
function αβγ_numeric(states_sol::AbstractSolution)
    s = structure(analysis(states_sol))
    # Node at (Lᵢ, Lⱼ, Lₖ)
    n₇ = nodes(s)[7]
    displacements_n₇ = displacements(states_sol_case₁, n₇)
    # Displacements in the x (component 1) axis at node 7
    numerical_uᵢ = displacements_n₇[1]
    numerical_α = 1 .+ numerical_uᵢ / Lᵢ
    # Displacements in the y (component 2) axis at node 7
    numerical_uⱼ = displacements_n₇[2]
    numerical_β = 1 .+ numerical_uⱼ / Lⱼ
    # Displacements in the z (component 3) axis at node 7
    numerical_uₖ = displacements_n₇[3]
    numerical_γ = 1 .+ numerical_uₖ / Lₖ
    return numerical_α, numerical_β, numerical_γ, numerical_uᵢ, numerical_uⱼ, numerical_uₖ
end
# Numeric solution for testing
numeric_α_case₁, numeric_β_case₁, numeric_γ_case₁, numeric_uᵢ_case₁, _, _ = αβγ_numeric(states_sol_case₁)
# Extract ℙ and ℂ from the last state using a random element
e = rand(elements(s₁))
# Cosserat or second Piola-Kirchhoff stress tensor
ℙ_numeric_case₁ = stress(states_sol_case₁, e)
# ℙᵢᵢ component: 
ℙᵢᵢ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 1, 1)
# ℙⱼⱼ component: 
ℙⱼⱼ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 2, 2)
# ℙₖₖ component: 
ℙₖₖ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 3, 3)
# Get the Right hand Cauchy strain tensor ℂ at a random state 
ℂ_rand_numeric_case₁ = rand(strain(states_sol_case₁, e))
# Get the Second Piola Kirchhoff stress tensor ℙ at a random state 
ℙ_rand_numeric_case₁ = rand(stress(states_sol_case₁, e))
# Load factors 
load_factors_case₁ = load_factors(sa₁)
# -----------------------------------------------
# Case 2 - GMSH mesh and `HyperElastic` material
#------------------------------------------------
# -------------------------------
# Materials
# -------------------------------
# Define a new HyperElastic material from the strain energy function
"Neo-Hookean strain energy function given the Green-Lagrange strain
tensor `𝔼`, second lamé parameter `μ` and bulk modulus `K`."
function strain_energy_neo(𝔼::AbstractMatrix, K::Real, μ::Real)
    # Right hand Cauchy strain tensor
    ℂ = Symmetric(2 * 𝔼 + eye(3))
    J = sqrt(det(ℂ))
    # First invariant
    I₁ = tr(ℂ)
    # Strain energy function 
    Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
end
params = [K, μ] # The order must be the same defined in the strain energy (splatting)
neo_hookean_hyper = HyperElastic(params, strain_energy_neo, "neoHyper")
# Material types without assigned elements
mat_types = [neo_hookean_hyper]
s_materials = StructuralMaterials(mat_types)
# -------------------------------
# Boundary Conditions
# -------------------------------
# Redefine the load boundary condition 
bc₄ = LocalPressureBoundaryCondition([:u], t -> [p * t], "tension")
# BoundaryConditions types without assigned node, feces and elements
vbc = [bc₁, bc₂, bc₃, bc₄]
s_boundary_conditions = StructuralBoundaryConditions(vbc)
# -------------------------------
# Entities
# -------------------------------
# Entities types without assigned nodes, faces and elements
vfaces = [TriangularFace("triangle")]
velems = [Tetrahedron("tetrahedron")]
s_entities = StructuralEntities(velems, vfaces)
# -------------------------------
# Mesh
# -------------------------------
file_name_msh = joinpath(@__DIR__, "uniaxial_compression.msh")
if GENERATE_MSH
    file_name_geo = joinpath(@__DIR__, "uniaxial_compression.geo")
    run(`gmsh -3 $file_name_geo -o $file_name`)
end
msh_file = MshFile(file_name_msh)
# -------------------------------
# Structure
# -------------------------------
s₂ = Structure(msh_file, s_materials, s_boundary_conditions, s_entities)
# Final load factor
sa₂ = StaticAnalysis(s₂, NSTEPS=NSTEPS)
reset!(sa₂)
# -------------------------------
# Numerical solution
# -------------------------------
# Extract ℙ and ℂ from the last state using a random element
e = rand(elements(s₂))
states_sol_case₂ = solve(sa₂, nr)
# Numeric solution for testing
numeric_α_case₂, numeric_β_case₂, numeric_γ_case₂, numeric_uᵢ_case₂, _, _ = αβγ_numeric(states_sol_case₂)
# Cosserat or second Piola-Kirchhoff stress tensor
ℙ_numeric_case₂ = stress(states_sol_case₂, e)
# ℙᵢᵢ component: 
ℙᵢᵢ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 1, 1)
# ℙⱼⱼ component: 
ℙⱼⱼ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 2, 2)
# ℙₖₖ component: 
ℙₖₖ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 3, 3)
# Get the Right hand Cauchy strain tensor ℂ at a random state 
ℂ_rand_numeric_case₂ = rand(strain(states_sol_case₂, e))
# Get the Second Piola Kirchhoff stress tensor ℙ at a random state 
ℙ_rand_numeric_case₂ = rand(stress(states_sol_case₂, e))
# Load factors 
load_factors_case₂ = load_factors(sa₂)
#-----------------------------
# Analytic solution  
#-----------------------------
# Test with Second Piola-Kirchoff stress tensor `ℙ`.
"Computes ℙ(1,1) given α, β and γ."
analytic_ℙᵢᵢ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K) =
    μ * α - μ * (α .^ (-1)) + K * (β .^ 2) .* (α .* (β .^ 2) .- 1)
"Computes ℙ(2,2) given α, β and γ."
analytic_ℙⱼⱼ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K) =
    μ * β - μ * (β .^ (-1)) + K * β .* ((α .^ 2) .* (β .^ 2) - α)
"Computes ℙ(2,2) given α, β and γ."
analytic_ℙₖₖ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K) =
    analytic_ℙⱼⱼ(α, β, μ, K)
# Compute the analytic Second Piola-Kirchoff stress tensor `ℙ` for the numeric vectors α and β
# Case 1 
ℙᵢᵢ_analytic_case₁ = analytic_ℙᵢᵢ(numeric_α_case₁, numeric_β_case₁)
ℙⱼⱼ_analytic_case₁ = analytic_ℙⱼⱼ(numeric_α_case₁, numeric_β_case₁)
ℙₖₖ_analytic_case₁ = analytic_ℙₖₖ(numeric_α_case₁, numeric_β_case₁)
# Case 2 
ℙᵢᵢ_analytic_case₂ = analytic_ℙᵢᵢ(numeric_α_case₂, numeric_β_case₂)
ℙⱼⱼ_analytic_case₂ = analytic_ℙⱼⱼ(numeric_α_case₂, numeric_β_case₂)
ℙₖₖ_analytic_case₂ = analytic_ℙₖₖ(numeric_α_case₂, numeric_β_case₂)
#-----------------------------
# Test boolean for CI  
#-----------------------------

@testset "Case 1 Uniaxial Compression Example" begin
    @test ℙᵢᵢ_analytic_case₁ ≈ ℙᵢᵢ_numeric_case₁ rtol = RTOL
    @test ℙᵢᵢ_analytic_case₁ ≈ ℙᵢᵢ_numeric_case₁ rtol = RTOL
    @test ℙⱼⱼ_analytic_case₁ ≈ ℙⱼⱼ_numeric_case₁ atol = ATOL
    @test ℙₖₖ_analytic_case₁ ≈ ℙₖₖ_numeric_case₁ atol = ATOL
    @test norm(ℙⱼⱼ_analytic_case₁) ≈ 0 atol = ATOL
    @test norm(ℙₖₖ_analytic_case₁) ≈ 0 atol = ATOL
    @test p * load_factors_case₁ ≈ ℙᵢᵢ_analytic_case₁ rtol = RTOL
end

@testset "Case 2 Uniaxial Compression Example" begin
    @test ℙᵢᵢ_analytic_case₂ ≈ ℙᵢᵢ_numeric_case₂ rtol = RTOL
    @test ℙᵢᵢ_analytic_case₂ ≈ ℙᵢᵢ_numeric_case₂ rtol = RTOL
    @test ℙⱼⱼ_analytic_case₂ ≈ ℙⱼⱼ_numeric_case₂ atol = ATOL
    @test ℙₖₖ_analytic_case₂ ≈ ℙₖₖ_numeric_case₂ atol = ATOL
    @test norm(ℙⱼⱼ_analytic_case₂) ≈ 0 atol = ATOL
    @test norm(ℙₖₖ_analytic_case₂) ≈ 0 atol = ATOL
    @test p * load_factors_case₂ ≈ ℙᵢᵢ_analytic_case₂ rtol = RTOL
end

