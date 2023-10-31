using Test, LinearAlgebra, SparseArrays
using Dictionaries: dictionary

using ONSAS.StaticAnalyses
using ONSAS.Circles
using ONSAS.SVKMaterial
using ONSAS.FixedDofBoundaryConditions
using ONSAS.GlobalLoadBoundaryConditions
using ONSAS.StructuralAnalyses
using ONSAS.Solvers
using ONSAS.Assemblers
using ONSAS.Solutions
using ONSAS.Structures
using ONSAS.StructuralBoundaryConditions
using ONSAS.StructuralMaterials
using ONSAS.Nodes
using ONSAS.Trusses
using ONSAS.StructuralSolvers
using ONSAS.Meshes
using ONSAS.StaticStates
using ONSAS.NonLinearStaticAnalyses

const RTOL = 5e-2

## scalar parameters
E = 210e9  # Young modulus in Pa
ν = 0.0  # Poisson's modulus
A = 2.5e-3   # Cross-section area in m^2
ang = 65 # truss angle in degrees
L = 2 # Length in m
d = L * cos(deg2rad(65))   # vertical distance in m
h = L * sin(deg2rad(65))
Fₖ = -3e8  # vertical load in N
Fⱼ = -2e8  # vertical load in N
# -------------------------------
# Materials
# -------------------------------
steel = SVK(; E=E, ν=ν, label="steel")
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
nodes = [n₁, n₂, n₃]
## Entities
truss₁ = Truss(n₁, n₂, s, "left_truss") # [n₁, n₂]
truss₂ = Truss(n₂, n₃, s, "right_truss") # [n₂, n₃]
elements = [truss₁, truss₂]
## Mesh
s_mesh = Mesh(; nodes, elements)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 3
set_dofs!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
mat_dict = dictionary([steel => [truss₁, truss₂]])
s_materials = StructuralMaterial(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁ = FixedDof(:u, [1, 2, 3], "fixed_uₓ_uⱼ_uₖ")
bc₂ = FixedDof(:u, [2], "fixed_uⱼ")
# Load
bc₃ = GlobalLoad(:u, t -> [0, 0, Fₖ * t], "load in k")
bc₄ = GlobalLoad(:u, t -> [0, Fⱼ * t, 0.0], "load in j")
node_bcs = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂], bc₄ => [n₂]])
s_boundary_conditions = StructuralBoundaryCondition(; node_bcs)
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
# Random static state
ΔUᵏ = rand(2)
Uᵏ = rand(9)
Fₑₓₜᵏ = rand(9)
Fᵢₙₜᵏ = rand(9)
Kₛᵏ = spzeros(9, 9)
ϵᵏ = dictionary([truss₁ => rand(3, 3), truss₂ => rand(3, 3)])
σᵏ = dictionary([truss₁ => rand(3, 3), truss₂ => rand(3, 3)])
s_assembler = Assembler(2)
iter_residuals = ResidualsIterationStep()
res_forces = zeros(2)

sst_rand = FullStaticState(free_dofs(s), ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ,
                           s_assembler,
                           iter_residuals)

@testset "ONSAS.StructuralAnalyses.StaticAnalyses.FullStaticState" begin

    # Accessors
    @test displacements(sst_rand) == Uᵏ
    @test Δ_displacements(sst_rand) == ΔUᵏ
    @test external_forces(sst_rand) == Fₑₓₜᵏ
    @test iteration_residuals(sst_rand) == iter_residuals
    @test residual_forces!(sst_rand) == Fₑₓₜᵏ[free_dofs(s)] - Fᵢₙₜᵏ[free_dofs(s)]
    @test tangent_matrix(sst_rand) == Kₛᵏ
    @test strain(sst_rand) == ϵᵏ
    @test stress(sst_rand) == σᵏ
    @test free_dofs(sst_rand) == free_dofs(s)

    # Iteration
    @test assembler(sst_rand) == s_assembler
    @test iteration_residuals(sst_rand) == iter_residuals
    norm_r = norm(residual_forces!(sst_rand))
    relative_norm_res = norm_r / norm(external_forces(sst_rand))
    @test residual_forces_norms(sst_rand) == (norm_r, relative_norm_res)
    norm_ΔU = norm(Δ_displacements(sst_rand))
    norm_U = norm(displacements(sst_rand))
    @test residual_displacements_norms(sst_rand) == (norm_ΔU, norm_ΔU / norm_U)

    ΔUᵏ⁺¹ = rand(2)
    sst_rand.Uᵏ[free_dofs(s)] += ΔUᵏ⁺¹
    Uᵏ⁺¹ = Uᵏ
    Uᵏ⁺¹[free_dofs(s)] += ΔUᵏ⁺¹
    @test displacements(sst_rand) == Uᵏ⁺¹

    # Reset the assembled magnitudes
    reset!(sst_rand)
    @test isempty(assembler(sst_rand).V)
    @test isempty(assembler(sst_rand).I)
    @test isempty(assembler(sst_rand).J)
    @test iszero(internal_forces(sst_rand))
    @test iszero(tangent_matrix(sst_rand))
    @test iszero(displacements(sst_rand))
    @test iszero(Δ_displacements(sst_rand))
    @test all([iszero(strain(sst_rand)[e]) for e in Meshes.elements(s)])
    @test all([iszero(stress(sst_rand)[e]) for e in Meshes.elements(s)])

    # Default static analysis of the structure
    default_s = FullStaticState(s)

    # Assemble process
    # truss₁ element
    fᵢₙₜ_e_1 = rand(6)
    k_e_1 = rand(6, 6)
    σ_e_1 = rand(3, 3)
    ϵ_e_1 = Symmetric(rand(3, 3))
    assemble!(default_s, fᵢₙₜ_e_1, truss₁)
    @test internal_forces(default_s)[1:6] ≈ fᵢₙₜ_e_1 rtol = RTOL
    assemble!(default_s, k_e_1, truss₁)
    @test internal_forces(default_s)[1:6] ≈ fᵢₙₜ_e_1 rtol = RTOL
    assemble!(default_s, σ_e_1, ϵ_e_1, truss₁)
    # truss₂ element
    fᵢₙₜ_e_2 = rand(6)
    k_e_2 = rand(6, 6)
    σ_e_2 = rand(3, 3)
    ϵ_e_2 = Symmetric(rand(3, 3))
    assemble!(default_s, fᵢₙₜ_e_2, truss₂)
    assemble!(default_s, k_e_2, truss₂)
    assemble!(default_s, σ_e_2, ϵ_e_2, truss₂)
    # End assemble
    end_assemble!(default_s)

    # Manufactured assemble
    Fᵢₙₜ = zeros(9)
    Fᵢₙₜ[1:6] += fᵢₙₜ_e_1
    Fᵢₙₜ[4:9] += fᵢₙₜ_e_2
    K_system = zeros(9, 9)
    K_system[1:6, 1:6] += k_e_1
    K_system[4:9, 4:9] += k_e_2

    @test internal_forces(default_s) ≈ Fᵢₙₜ rtol = RTOL
    @test tangent_matrix(default_s) ≈ K_system rtol = RTOL
    @test strain(default_s)[truss₁] == ϵ_e_1
    @test strain(default_s)[truss₂] == ϵ_e_2
    @test stress(default_s)[truss₁] == σ_e_1
    @test stress(default_s)[truss₂] == σ_e_2
end

# NonLinearStaticAnalysis with a final load factor
λ₁ = 10
NSTEPS = 9
init_step = 7

sa = NonLinearStaticAnalysis(s, λ₁; NSTEPS=NSTEPS)
sa_init = NonLinearStaticAnalysis(s, λ₁; NSTEPS=NSTEPS, initial_step=init_step)

@testset "ONSAS.StructuralAnalyses.StaticAnalyses.NonLinearStaticAnalysis" begin
    @test structure(sa_init) == s
    λ₀ = λ₁ / NSTEPS
    λᵥ = LinRange(λ₀, λ₁, NSTEPS)
    @test initial_time(sa_init) == first(λᵥ)
    @test current_time(sa_init) == λᵥ[init_step]
    @test final_time(sa_init) == last(λᵥ)
    @test current_load_factor(sa_init) == init_step * λ₁ / NSTEPS
    @test load_factors(sa_init) == λᵥ

    # External forces
    @test external_forces(current_state(sa_init)) == zeros(num_dofs(s))
    apply!(sa_init, load_bcs(s_boundary_conditions))
    @test external_forces(current_state(sa_init))[4:6] ==
          bc₄(current_time(sa_init)) + bc₃(current_time(sa_init))

    # Next step
    next!(sa_init)
    @test current_time(sa_init) == λᵥ[init_step + 1]
    @test current_load_factor(sa_init) == (init_step + 1) * λ₁ / NSTEPS
    @test !is_done(sa_init)
    next!(sa_init)
    next!(sa_init)
    @test is_done(sa_init)

    # Reset the analysis
    reset!(sa_init)
    @test current_time(sa_init) == first(λᵥ)
end

@testset "ONSAS.StructuralSolvers.Solution" begin
    solved_states = [sst_rand, sst_rand, sst_rand]
    num_states = length(solved_states)
    nr = NewtonRaphson()
    states_sol = Solution(FullStaticState[], sa, nr)
    foreach(solved_states) do st
        push!(states_sol, st)
    end

    @test length(states(states_sol)) == num_states
    @test analysis(states_sol) == sa
    @test solver(states_sol) == nr

    dof = Dof(5)
    vdof = [Dof(1), Dof(4)]

    @test displacements(states_sol, dof) == repeat([displacements(sst_rand)[dof]], num_states)

    @test displacements(states_sol, vdof) ==
          [repeat([displacements(sst_rand)[vdof[1]]], num_states),
           repeat([displacements(sst_rand)[vdof[2]]], num_states)]

    @test displacements(states_sol, n₁) ==
          [repeat([displacements(sst_rand)[dofs(n₁)[:u][1]]], num_states),
           repeat([displacements(sst_rand)[dofs(n₁)[:u][2]]], num_states),
           repeat([displacements(sst_rand)[dofs(n₁)[:u][3]]], num_states)]

    @test internal_forces(states_sol, dof) == repeat([internal_forces(sst_rand)[dof]], num_states)

    @test internal_forces(states_sol, vdof) ==
          [repeat([internal_forces(sst_rand)[vdof[1]]], num_states),
           repeat([internal_forces(sst_rand)[vdof[2]]], num_states)]

    @test internal_forces(states_sol, n₁) ==
          [repeat([internal_forces(sst_rand)[dofs(n₁)[:u][1]]], num_states),
           repeat([internal_forces(sst_rand)[dofs(n₁)[:u][2]]], num_states),
           repeat([internal_forces(sst_rand)[dofs(n₁)[:u][3]]], num_states)]

    @test external_forces(states_sol, dof) == repeat([external_forces(sst_rand)[dof]], num_states)

    @test external_forces(states_sol, vdof) ==
          [repeat([external_forces(sst_rand)[vdof[1]]], num_states),
           repeat([external_forces(sst_rand)[vdof[2]]], num_states)]

    @test external_forces(states_sol, n₁) ==
          [repeat([external_forces(sst_rand)[dofs(n₁)[:u][1]]], num_states),
           repeat([external_forces(sst_rand)[dofs(n₁)[:u][2]]], num_states),
           repeat([external_forces(sst_rand)[dofs(n₁)[:u][3]]], num_states)]

    iteration_residuals(states_sol)
end
