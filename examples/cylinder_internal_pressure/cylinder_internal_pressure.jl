# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example  
#----------------------------------------------------
using LinearAlgebra: norm
using ONSAS
using StaticArrays: SVector
using Test: @test, @testset
using Suppressor: @capture_out
using TimerOutputs: @timeit, reset_timer!, print_timer
## scalar parameters (dimensions in mm an MPa)
Lₖ = 30; # cylinder length in 𝐞ₖ mm
Rᵢ = 100; # inner radius in mm
Rₑ = 200; # outer radius in mm
p = 30; # internal pressure in MPa
E = 210.0;  # Young modulus in MPa
ν = 0.3;  # Poisson ratio
pressure(t::Real) = -p * t
## number of steps 
NSTEPS = 9
## tolerances for testing
ATOL = 1e-2 * (Rₑ - Rᵢ)
## Plot results
plot_results = true
"Creates the mesh."
function create_msh()
    # Run gmsh to generate the mesh
    command = `gmsh -3 examples/cylinder_internal_pressure/cylinder.geo`
    # Generate the mesh and capture the outputs
    o = @capture_out begin
        run(command)
    end
    gmsh_println(o)
end
out_gmsh = create_msh()
"Builds the `Structure`."
function cylinder_structure(
    material::AbstractMaterial,
    L::Real, Rᵢ::Real, Rₑ::Real,
    p::Function
)
    # Materials
    # -------------------------------
    materials = StructuralMaterials(material)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    # Dirichlet boundary conditions 
    bc₁_label = "fixed-ui"
    bc₂_label = "fixed-uj"
    bc₃_label = "fixed-uk"
    bc₄_label = "pressure"
    bc₁ = FixedDofBoundaryCondition([:u], [1], bc₁_label)
    bc₂ = FixedDofBoundaryCondition([:u], [2], bc₂_label)
    bc₃ = FixedDofBoundaryCondition([:u], [3], bc₃_label)
    # Neumann boundary conditions 
    bc₄ = LocalPressureBoundaryCondition([:u], t -> p(t), bc₄_label)
    boundary_conditions = StructuralBoundaryConditions(bc₁, bc₂, bc₃, bc₄)
    bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
    # -------------------------------
    # Entities
    # -------------------------------
    # Entities types without assigned nodes, faces and elements
    faces_label = "triangle"
    elements_label = "tetrahedron"
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elements_label)]
    entities = StructuralEntities(velems, vfaces)
    entities_labels = [faces_label, elements_label]
    # -------------------------------
    # Mesh
    # -------------------------------
    msh_path = joinpath(@__DIR__, "cylinder.msh")
    msh_mesh = MshFile(msh_path)
    # -------------------------------
    # Structure
    # -------------------------------
    structure = Structure(msh_mesh, materials, boundary_conditions, entities)
end;
# -------------------------------
# Structure
# -------------------------------
reset_timer!()
# -------------------------------
# Materials
# -------------------------------
mat_label = "mat";
linear_material = IsotropicLinearElastic(E, ν, mat_label);
@timeit "Building the linear structure ⚪ " begin
    linear_cylinder = cylinder_structure(linear_material, Lₖ, Rᵢ, Rₑ, pressure)
end
svk_material = SVK(E=E, ν=ν, label=mat_label);
@timeit "Building the non-linear structure 🔘" begin
    nonlinear_cylinder = cylinder_structure(svk_material, Lₖ, Rᵢ, Rₑ, pressure)
end
# -------------------------------
# Structural Analysis
# -------------------------------
"Defines an structural `AbstractStaticAnalysis`."
static_analysis(structure::Structure, analysis::Type{<:AbstractStaticAnalysis}; NSTEPS::Int) =
    analysis(structure, NSTEPS=NSTEPS);
# -----------------------------------------------
# Case 1 - Static linear elastic case
#------------------------------------------------
@timeit "Defining the linear analysis 👷 🔎 ⚪" begin
    linear_analysis = static_analysis(linear_cylinder, LinearStaticAnalysis, NSTEPS=NSTEPS)
end
# -----------------------------------------------
# Case 2 - Static non-linear elastic case
#----------------------------------------------
@timeit "Defining the non-linear analysis 👲 🔎 🔘" begin
    nonlinear_analysis = static_analysis(nonlinear_cylinder, NonLinearStaticAnalysis, NSTEPS=NSTEPS)
end
# -------------------------------
# Numerical solution
# -------------------------------
# Linear analysis
# -------------------------------
@timeit "Solving the linear analysis 🐢->🐕" begin
    states_lin_sol = solve!(linear_analysis)
end
# get time vector or load factors  
λᵥ = load_factors(linear_analysis)
# Get the solution at a random point 
"Return a rand point in the cylinder (R, θ, L)."
rand_point_cylinder(Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ, Lₖ::Real=Lₖ) =
    [rand() * (Rₑ - Rᵢ) + Rᵢ, rand() * 2 * π, rand() * Lₖ]
# Get the internal radial displacement at p = (0, Rᵢ, 0)
linear_cylinder_mesh = mesh(linear_cylinder)
nᵢ = nodes(linear_cylinder_mesh)[4];
uᵣ_numeric_nᵢ = displacements(states_lin_sol, nᵢ, 2);
# Get the external radial displacement at p = (-Rₑ, 0, Lₖ)
nₑ = nodes(linear_cylinder_mesh)[15];
uᵣ_numeric_nₑ = displacements(states_lin_sol, nₑ, 1);
# Generate a random point
rand_R, rand_θ, rand_z = rand_point_cylinder();
# rand_R = Rᵢ;
# rand_θ = 0.0
p_rand = SVector(rand_R * cos(rand_θ), rand_R * sin(rand_θ), rand_z);
# Displacements at p
point_evaluator = PointEvalHandler(linear_cylinder_mesh, p_rand);
uᵢ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 1);
uⱼ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 2);
uₖ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 3);
# 
uᵣ_numeric_p_rand = sqrt.(@. uᵢ_numeric_p_rand^2 + uⱼ_numeric_p_rand^2);
#  Non-linear analysis
# -------------------------------
tols = ConvergenceSettings(rel_U_tol=1e-8, rel_res_force_tol=1e-8, max_iter=30)
alg = NewtonRaphson(tols)
@timeit "Solving the non-linear analysis 🐢->🐕" begin
    states_nonlinear_sol = solve!(nonlinear_analysis, alg)
end
# Get the internal radial displacement at p = (0, Rᵢ, 0)
uᵣ_numeric_nonlinear_nᵢ = displacements(states_nonlinear_sol, nᵢ, 2);
# Get the external radial displacement at p = (-Rₑ, 0, Lₖ)
uᵣ_numeric_nonlinear_nₑ = displacements(states_nonlinear_sol, nₑ, 1);
# -------------------------------
# Analytic solution
# -------------------------------
t = last(λᵥ)
"Analytic radial displacements uᵣ at radius`r` and time `t`."
function uᵣ(
    r::Real, t::Real,
    E::Real=E, ν::Real=ν, p::Function=pressure,
    Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ,
)
    "Constant A for the analytic solution."
    A(t::Real, Rᵢ::Real, Rₑ::Real, E::Real, ν::Real, p::Function) =
        (1 + ν) * (1 - 2 * ν) * Rᵢ^2 * -p(t) / (E * (Rₑ^2 - Rᵢ^2))
    "Constant B for the analytic solution."
    B(t::Real, Rᵢ::Real, Rₑ::Real, E::Real, ν::Real, p::Function) =
        (1 + ν) * Rᵢ^2 * Rₑ^2 * -p(t) / (E * (Rₑ^2 - Rᵢ^2))
    uᵣ = A(t, Rᵢ, Rₑ, E, ν, p) * r + B(t, Rᵢ, Rₑ, E, ν, p) / r
end;
uᵣ_analytic_nᵢ = [uᵣ(Rᵢ, t) for t in λᵥ];
uᵣ_analytic_nₑ = [uᵣ(Rₑ, t) for t in λᵥ];
uᵣ_analytic_p_rand = [uᵣ(rand_R, t) for t in λᵥ];
#-----------------------------
# Print & plots 
#-----------------------------
print_timer(title="Analysis with $(out_gmsh[:nelems]) elements and $(out_gmsh[:nnodes]) nodes")
if plot_results
    using Plots
    vec_p = [-pressure(λ) for λ in λᵥ]
    vec_p_non_in = [-pressure(λ) for λ in load_factors(nonlinear_analysis)]
    fig = plot(
        vec_p, uᵣ_numeric_nᵢ, label="numeric linear uᵣ n=(0, Rᵢ, 0)",
        legend=:topleft, color=:orange, lw=2, ls=:dash, markershape=:circle,
    )
    plot!(fig,
        vec_p, -uᵣ_numeric_nₑ, label="numeric linear uᵣ n=(-Rₑ, 0 , Lₖ)",
        legend=:topleft, color=:skyblue, lw=2, ls=:solid, markershape=:square,
    )
    plot!(fig,
        vec_p, uᵣ_analytic_nᵢ, label="analytic linear uᵣ(Rᵢ)",
        legend=:topleft, color=:black, lw=2, ls=:dash, markershape=:none,
    )
    plot!(fig,
        vec_p, uᵣ_analytic_nₑ, label="analytic linear uᵣ(Rₑ)",
        legend=:topleft, color=:black, lw=2, ls=:solid
    )
    # Plot comparing linear and non linear solutions 
    plot!(fig,
        vec_p_non_in, uᵣ_numeric_nonlinear_nᵢ, label="non-linear uᵣ(0, Rᵢ, 0)",
        color=:red, lw=2, marker=:circle, markersize=3
    )
    plot!(fig,
        vec_p_non_in, -uᵣ_numeric_nonlinear_nₑ, label="non-linear uᵣ(-Rₑ, 0 , Lₖ)",
        color=:blue, lw=2, marker=:circle, markersize=3
    )
    # add labels
    xlabel!("λᵥ [MPa]")
    ylabel!("uᵣ [mm]")
    display(fig)
end
#-----------------------------
# Test booleans
#-----------------------------
# Test symmetry and boundary conditions for a random slice
#-------------------------------------------
function test_solution_at_slice(
    sol::AbstractSolution=states_lin_sol;
    atol::Real=ATOL, atolr=ATOLR,
    Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ, Lₖ::Real=Lₖ
)
    structure = ONSAS.structure(analysis(sol))
    # Generic surface s at z = Lₖ 
    rand_R, rand_θ₁, Lₖ = rand_point_cylinder(Rᵢ, Rₑ, Lₖ)
    # Set by force Lₖ
    rand_θ₂ = rand() * 2 * π
    # Random point ∈ axis x
    p_randᵢ = [rand_R, 0.0, Lₖ]
    # Random point ∈ axis y
    p_randⱼ = [0.0, rand_R, Lₖ]
    # Random point between the internal and external surface
    p_rand₁ = [rand_R * cos(rand_θ₁), rand_R * sin(rand_θ₁), Lₖ]
    p_rand₂ = [rand_R * cos(rand_θ₂), rand_R * sin(rand_θ₂), Lₖ]
    # Vector of points to test
    vec_points = [p_randᵢ, p_randⱼ, p_rand₁, p_rand₂]
    #
    point_evaluator = PointEvalHandler(mesh(structure), vec_points)
    U = displacements(sol, point_evaluator)
    # Check uₖ = 0 ∀ p ∈ s
    zero_uₖ = all([≈(norm(u[3]), 0.0, atol=atol) for u in U])
    # Check uᵢ = 0 ∀ p ∈ s & ∈ axis y
    Uᵢ_in_axis_y = displacements(states_lin_sol, point_evaluator, 2)[1]
    zero_uₖ_axis_y = all([≈(norm(uᵢ_p_in_axis_y), 0.0, atol=atol) for uᵢ_p_in_axis_y in Uᵢ_in_axis_y])
    # Check uⱼ = 0 ∀ p ∈ s & ∈ axis x
    Uⱼ_in_axis_x = displacements(states_lin_sol, point_evaluator, 1)[2]
    zero_uⱼ_axis_x = all([≈(norm(uⱼ_p_in_axis_y), 0.0, atol=atol) for uⱼ_p_in_axis_y in Uⱼ_in_axis_x])
    # Check uᵣ(r,θ₁) =  uᵣ(r,θ₁)  at last time
    rand₁_index = 3
    uᵣ_rand₁ = sum(last.(U[rand₁_index][1:2]) .^ 2)
    rand₂_index = 4
    uᵣ_rand₂ = sum(last.(U[rand₂_index][1:2]) .^ 2)
    uᵣ_not_depends_on_θ = ≈(uᵣ_rand₁, uᵣ_rand₂, atol=atolr)
    return uᵣ_not_depends_on_θ, zero_uₖ, zero_uₖ_axis_y, zero_uⱼ_axis_x
end;
# Test simmetry and boundary conditions
uᵣ_not_depends_on_θ, zero_uₖ, zero_uₖ_axis_y, zero_uⱼ_axis_x =
    test_solution_at_slice(states_lin_sol, atol=ATOL, atolr=10 * ATOL)
@testset "Case 1: Linear Analysis " begin
    @info "uᵣ(r,θ₁,L₁) = uᵣ(r,θ₂,L₂)?" uᵣ_not_depends_on_θ
    @test zero_uₖ
    @test zero_uₖ_axis_y
    @test zero_uⱼ_axis_x
    @test uᵣ_numeric_p_rand ≈ uᵣ_analytic_p_rand atol = ATOL
    @test uᵣ_analytic_nᵢ ≈ uᵣ_numeric_nᵢ atol = ATOL
    @test uᵣ_analytic_nₑ ≈ -uᵣ_numeric_nₑ atol = ATOL
end
# Test simmetry and boundary conditions
uᵣ_not_depends_on_θ_case2, zero_uₖ_case2, zero_uₖ_axis_y_case2, zero_uⱼ_axis_x_case2 =
    test_solution_at_slice(states_lin_sol, atol=ATOL, atolr=10 * ATOL)
@testset "Case 2: Non-Linear Analysis " begin
    @info "uᵣ(r,θ₁,L₁) = uᵣ(r,θ₂,L₂)?" uᵣ_not_depends_on_θ_case2
    @test zero_uₖ_case2
    @test zero_uₖ_axis_y_case2
    @test zero_uⱼ_axis_x_case2
end
