
using Test: @testset, @test

using ONSAS

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
add_dofs!(s_mesh, :u, 3)
add_dofs!(s_mesh, :θ, 3)
# -------------------------------
# Materials
# -------------------------------
mat_dict = dictionary([steel => [truss₁, truss₂]])
s_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
bc₁ = PinnedDisplacementBoundaryCondition("pinned")
bc₂ = FⱼLoadBoundaryCondition(Fⱼ, "load in j")
node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂]])
s_boundary_conditions = StructuralBoundaryConditions(node_bc)
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


@testset "ONSAS.StructuralAnalyses.StaticState" begin

    # Random static state
    Δuᵏ = zeros(2)
    Uᵏ = zeros(4)
    Fₑₓₜᵏ = rand(4)
    Fᵢₙₜᵏ = rand(4)
    Kₛᵏ = rand(4, 4)
    ϵ = Vector{}(undef, 2)
    σ = Vector{}(undef, 2)
    assembler = Assembler(2)

    sa_rand = StaticState(Δuᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵ, σ, assembler)
    @test residual_forces(sa_rand) == Fₑₓₜᵏ - Fᵢₙₜᵏ
    @test tangent_matrix(sa_rand) == Kₛᵏ

end


@testset "ONSAS.StructuralAnalyses.StaticAnalysis" begin
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