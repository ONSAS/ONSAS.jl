using ONSAS.StructuralAnalyses.StaticAnalyses
using Test

## scalar parameters
E = 210e9  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 2.5e-3   # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m 
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
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
add!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
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
s_boundary_conditions = StructuralBoundaryConditions(node_bc)
# Simplificar en un solo diccionario y que el constructor reciba las dos
# -------------------------------
# Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_boundary_conditions)
# -------------------------------
# Solver
# -------------------------------
tol_f = 1e-10;
tol_u = 1e-10;
max_iter = 100;
tols = ConvergenceSettings(tol_f, tol_u, max_iter)
nr = NewtonRaphson(tols)

@testset "ONSAS.StructuralAnalyses.StaticAnalyses.StaticState" begin

    # Random static state
    Δuᵏ = zeros(2)
    Uᵏ = zeros(4)
    Fₑₓₜᵏ = rand(4)
    Fᵢₙₜᵏ = rand(4)
    Kₛᵏ = rand(4, 4)
    ϵ = Vector{}(undef, 2)
    σ = Vector{}(undef, 2)
    s_assembler = Assembler(2)
    iter_residuals = ResidualsIterationStep()

    sst_rand = StaticState(s, Δuᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵ, σ,
        s_assembler, iter_residuals)

    @test assembler(sst_rand) == s_assembler
    @test displacements(sst_rand) == Uᵏ
    @test Δ_displacements(sst_rand) == Δuᵏ
    @test external_forces(sst_rand) == Fₑₓₜᵏ
    @test internal_forces(sst_rand) == iter_residuals
    @test residual_forces(sst_rand) == Fₑₓₜᵏ[free_dofs(s)] - Fᵢₙₜᵏ[free_dofs(s)]
    @test tangent_matrix(sst_rand) == Kₛᵏ

end
#=

@testset "ONSAS.StructuralAnalyses..StaticAnalyses.StaticAnalysis" begin
    # Final load factor
    λ₁ = 10
    NSTEPS = 9
    init_step = 7
    sa_init = StaticAnalysis(s, λ₁, NSTEPS=NSTEPS, initial_step=init_step)

    @test structure(sa_init) == s
    λ₀ = λ₁ / NSTEPS
    λᵥ = LinRange(λ₀, λ₁, NSTEPS)
    @test initial_time(sa_init) == first(λᵥ)
    @test current_time(sa_init) == λᵥ[init_step]
    @test final_time(sa_init) == last(λᵥ)
    @test current_load_factor(sa_init) == init_step * λ₁ / NSTEPS
    @test load_factors(sa_init) == λᵥ

    # Next step 
    next!(sa_init)
    @test current_time(sa_init) == λᵥ[init_step+1]
    @test current_load_factor(sa_init) == (init_step + 1) * λ₁ / NSTEPS
    @test !is_done(sa_init)
    next!(sa_init)
    @test is_done(sa_init)

    # Reset and solve 
    sa = StaticAnalysis(s, NSTEPS=NSTEPS)

    solve(sa, nr)
end

=#