# ---------------------------------------------------------------- 
# Uniaxial Extension Example 1  from (Zerpa et. Al., 2019, CMAME).
# ----------------------------------------------------------------
using ONSAS.StaticAnalyses
using Test: @test, @testset
using StaticArrays: SVector
using LinearAlgebra: norm
## scalar parameters
E = 2.0                          # Young modulus in Pa
ν = 0.4                          # Poisson's ratio
p = 3                            # Tension load in Pa
const tension(t) = p * t         # Tension load function
Lᵢ = 2.0                         # Dimension in x of the box in m 
Lⱼ = 1.0                         # Dimension in y of the box in m
Lₖ = 1.0                         # Dimension in z of the box in m
const RTOL = 1e-4                # Relative tolerance for tests
const NSTEPS = 9                 # Number of steps for the test
include("./linear_extension_mesh.jl")   # Mesh Cube with Gmsh.jl
# -------------------------------
# Materials
# -------------------------------
mat_label = "mat"
mat = IsotropicLinearElastic(E, ν, mat_label)
# mat = SVK(E=E, ν=ν, label=mat_label)
s_materials = StructuralMaterials([mat])
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁_label = "fixed-ux"
bc₁ = FixedDofBoundaryCondition([:u], [1], bc₁_label)
bc₂_label = "fixed-uj"
bc₂ = FixedDofBoundaryCondition([:u], [2], bc₂_label)
bc₃_label = "fixed-uk"
bc₃ = FixedDofBoundaryCondition([:u], [3], bc₃_label)
# Load
bc₄_label = "tension"
bc₄ = GlobalLoadBoundaryCondition([:u], t -> [tension(t), 0, 0], bc₄_label)
# Create the  bcs vector
vbc = [bc₁, bc₂, bc₃, bc₄]
bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
s_boundary_conditions = StructuralBoundaryConditions(vbc)
# -------------------------------
# Entities
# -------------------------------
# Entities types without assigned nodes, faces and elements
faces_label = "triangle"
elems_label = "tetrahedron"
vfaces = [TriangularFace(faces_label)]
velems = [Tetrahedron(elems_label)]
entities_labels = [faces_label, elems_label]
s_entities = StructuralEntities(velems, vfaces)
# -------------------------------
# Mesh
# -------------------------------
filename = "linear_extension"
labels = [mat_label, entities_labels, bc_labels]
file_name_mesh = create_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename)
msh_file = MshFile(file_name_mesh)
# -------------------------------
# Structure
# -------------------------------
s = Structure(msh_file, s_materials, s_boundary_conditions, s_entities)
# -------------------------------
# Structural Analysis
# -------------------------------
sa = LinearStaticAnalysis(s, NSTEPS=NSTEPS)
# -------------------------------
# Numerical solution
# -------------------------------
states_sol = solve(sa)
# Select random points to test the solution
## Displacements
x₀_rand = Lᵢ * rand(2)
y₀_rand = Lⱼ * rand(2)
z₀_rand = Lₖ * rand(2)
p₁ = SVector(x₀_rand[1], y₀_rand[1], z₀_rand[1])
p₂ = SVector(x₀_rand[2], y₀_rand[2], z₀_rand[2])
# Evaluate the solution at p₁, p₂
eval_handler_rand = PointEvalHandler(mesh(s), [p₁, p₂])
# rand points displacements
# point 1
uᵢ_numeric_p₁ = displacements(states_sol, eval_handler_rand, 1)[1]
uⱼ_numeric_p₁ = displacements(states_sol, eval_handler_rand, 2)[1]
uₖ_numeric_p₁ = displacements(states_sol, eval_handler_rand, 3)[1]
# point 2
uᵢ_numeric_p₂ = displacements(states_sol, eval_handler_rand, 1)[2]
uⱼ_numeric_p₂ = displacements(states_sol, eval_handler_rand, 2)[2]
uₖ_numeric_p₂ = displacements(states_sol, eval_handler_rand, 3)[2]
## Strain and stresses
# Evaluate the solution at a random element
e_rand = rand(elements(s))
ϵ_e_rand = stress(states_sol, e_rand)
ϵᵢ_numeric_e_rand = getindex.(ϵ_e_rand, 1)
ϵⱼ_numeric_e_rand = getindex.(ϵ_e_rand, 2)
ϵₖ_numeric_e_rand = getindex.(ϵ_e_rand, 2)
# -------------------------------
# Analytic solution
# -------------------------------
## Displacements
"Computes displacements numeric solution uᵢ, uⱼ and uₖ for analytic validation."
function u_ijk_analytic(λᵥ::Vector{<:Real}, x₀::Real, y₀::Real, z₀::Real, ν::Real=ν, E::Real=E)

    𝐶(t) = tension(t) * (1 - ν - 2ν^2) / (1 - ν)

    uᵢ(t) = 𝐶(t) / E * x₀
    uⱼ(t) = 0.0
    uₖ(t) = 0.0

    return [[uᵢ(t) for t in λᵥ], [uⱼ(t) for t in λᵥ], [uₖ(t) for t in λᵥ]]
end
# point 1
u_analytic_p₁ = u_ijk_analytic(load_factors(sa), p₁[1], p₁[2], p₁[3])
uᵢ_analytic_p₁ = u_analytic_p₁[1]
uⱼ_analytic_p₁ = u_analytic_p₁[2]
uₖ_analytic_p₁ = u_analytic_p₁[3]
# point 2
u_analytic_p₂ = u_ijk_analytic(load_factors(sa), p₂[1], p₂[2], p₂[3])
uᵢ_analytic_p₂ = u_analytic_p₂[1]
uⱼ_analytic_p₂ = u_analytic_p₂[2]
uₖ_analytic_p₂ = u_analytic_p₂[3]
## Strains
"Computes strains numeric solution ϵᵢ, ϵⱼ and ϵₖ for analytic validation."
function ϵ_ijk_analytic(λᵥ::Vector{<:Real}, x₀::Real, y₀::Real, z₀::Real, ν::Real=ν, E::Real=E)

    𝐶(t) = tension(t) * (1 - ν - 2ν^2) / (1 - ν)

    ϵᵢ(t) = 𝐶(t) / E
    ϵⱼ(t) = 0.0
    ϵₖ(t) = 0.0

    return [[ϵᵢ(t) for t in λᵥ], [ϵⱼ(t) for t in λᵥ], [ϵₖ(t) for t in λᵥ]]
end
# point in the rand element selected
p_rand_e = rand(coordinates(e_rand))
ϵ_analytic_p_rand_e = ϵ_ijk_analytic(load_factors(sa), p_rand_e[1], p_rand_e[2], p_rand_e[3])
ϵᵢ_analytic_p_rand_e = ϵ_analytic_p_rand_e[1]
ϵⱼ_analytic_p_rand_e = ϵ_analytic_p_rand_e[2]
ϵₖ_analytic_p_rand_e = ϵ_analytic_p_rand_e[3]
#-----------------------------
# Test boolean for CI  
#-----------------------------
@testset "Linear Extension example" begin
    # Displacements
    @test uᵢ_numeric_p₁ ≈ uᵢ_analytic_p₁ rtol = RTOL
    @test norm(uⱼ_numeric_p₁) ≈ 0 atol = RTOL
    @test norm(uₖ_numeric_p₁) ≈ 0 atol = RTOL
    @test uᵢ_numeric_p₂ ≈ uᵢ_analytic_p₂ rtol = RTOL
    @test norm(uⱼ_numeric_p₂) ≈ 0 atol = RTOL
    @test norm(uₖ_numeric_p₂) ≈ 0 atol = RTOL
    # Strains
    @test ϵᵢ_numeric_e_rand ≈ ϵᵢ_e_rand rtol = RTOL
    @test norm(ϵⱼ_numeric_e_rand) ≈ 0 atol = RTOL
    @test norm(ϵₖ_numeric_e_rand) ≈ 0 atol = RTOL
end