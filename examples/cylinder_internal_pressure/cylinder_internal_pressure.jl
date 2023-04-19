# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example  
#----------------------------------------------------
using LinearAlgebra: norm
using ONSAS
using StaticArrays: SVector
using Test: @test, @testset
## scalar parameters (dimensions in mm an MPa)
L = 30; # cylinder length in 𝐞ₖ mm
Rᵢ = 100; # inner radius in mm
Rₑ = 200; # outer radius in mm
p = 10; # internal pressure in MPa
pressure(t::Real) = -p * t
E = 210.0;  # Young modulus in MPa
ν = 0.3;  # Poisson ratio


## tolerances for testing
ATOL = (Rₑ - Rᵢ) / 100

# Run gmsh to generate the mesh
command = `gmsh -3 examples/cylinder_internal_pressure/cylinder.geo`;
run(command);

"Builds the `Structure`."
function cylinder_structure(E::Real, ν::Real, L::Real, Rᵢ::Real, Rₑ::Real, p::Function)
    # -------------------------------
    # Materials
    # -------------------------------
    mat_label = "mat"
    material = IsotropicLinearElastic(E, ν, mat_label)
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
structure = cylinder_structure(E, ν, L, Rᵢ, Rₑ, pressure);
# -------------------------------
# Structural Analysis
# -------------------------------
"Defines an structural `AbstractStaticAnalysis`."
function static_analysis(
    structure::Structure, analysis::Type{<:AbstractStaticAnalysis};
    NSTEPS::Int
)
    analysis(structure, NSTEPS=NSTEPS)
end;
linear_analysis = static_analysis(structure, LinearStaticAnalysis, NSTEPS=2);
# -------------------------------
# Numerical solution
# -------------------------------
states_lin_sol, time, _ = @timed solve!(linear_analysis);
# 
λᵥ = load_factors(linear_analysis)
# Select a face to test the boundary conditions 
"Tests boundary conditions and symmetry for a slice at x = Lᵢ."
function test_bcs_slice(Lᵢ::Real; atol::Real=ATOL)
    # Generic surface s at x = Lᵢ 
    rand_R = rand() * (Rₑ - Rᵢ) + Rᵢ
    rand_θ₁ = rand() * 2 * π
    rand_θ₂ = rand() * 2 * π
    # Random point ∈ axis x
    p_randᵢ = [rand_R, 0.0, Lᵢ]
    # Random point ∈ axis y
    p_randⱼ = [0.0, rand_R, Lᵢ]
    # Random point between the internal and external surface
    p_rand₁ = [rand_R * cos(rand_θ₁), rand_R * sin(rand_θ₁), Lᵢ]
    p_rand₂ = [rand_R * cos(rand_θ₂), rand_R * sin(rand_θ₂), Lᵢ]
    # Vector of points to test
    vec_points = [p_randᵢ, p_randⱼ, p_rand₁, p_rand₂]
    #
    point_evaluator = PointEvalHandler(mesh(structure), vec_points)
    U = displacements(states_lin_sol, point_evaluator)
    # Check uₖ = 0 ∀ p ∈ s
    @test all([≈(norm(u[3]), 0.0, atol=atol) for u in U])
    # Check uᵢ = 0 ∀ p ∈ s & ∈ axis y
    Uᵢ_in_axis_y = displacements(states_lin_sol, point_evaluator, 2)[1]
    @test all([≈(norm(uᵢ_p_in_axis_y), 0.0, atol=atol) for uᵢ_p_in_axis_y in Uᵢ_in_axis_y])
    # Check uⱼ = 0 ∀ p ∈ s & ∈ axis x
    Uⱼ_in_axis_x = displacements(states_lin_sol, point_evaluator, 1)[2]
    @test all([≈(norm(uⱼ_p_in_axis_y), 0.0, atol=atol) for uⱼ_p_in_axis_y in Uⱼ_in_axis_x])
    # Check uᵣ(r,θ₁) =  uᵣ(r,θ₁)  at last time
    rand₁_index = 3
    @show uᵣ_rand₁ = sum(last.(U[rand₁_index][1:2]) .^ 2)
    rand₂_index = 4
    @show uᵣ_rand₂ = sum(last.(U[rand₂_index][1:2]) .^ 2)
end;
# Rand slice z coordiante
Lₖ = L / (1 + rand());
test_bcs_slice(Lₖ, atol=ATOL)

@info "The execution time was $(time) seconds 🔧."
# Displacements at internal and external nodes 
nᵢ = nodes(mesh(structure))[1]
@show uᵣᵢ_numeric = last(displacements(states_lin_sol, nᵢ, 1))
nₑ = nodes(mesh(structure))[8]
@show uᵣₑ_numeric = last(displacements(states_lin_sol, nₑ, 2))
# Displacements at a random point 
rand_R = rand() * (Rₑ - Rᵢ) + Rᵢ;
rand_θ = rand() * 2 * π;
p_rand = SVector(rand_R * cos(rand_θ), rand_R * sin(rand_θ), Lₖ);
point_evaluator = PointEvalHandler(mesh(structure), p_rand);
uᵢ_p_rand = displacements(states_lin_sol, point_evaluator, 1);
uⱼ_p_rand = displacements(states_lin_sol, point_evaluator, 2);
uₖ_p_rand = displacements(states_lin_sol, point_evaluator, 3);
uᵣ_p_rand = @. uᵢ_p_rand^2 + uⱼ_p_rand^2
# -------------------------------
# Analytic solution
# -------------------------------
t = last(λᵥ)
A(t, Rᵢ, Rₑ, E, ν) = (1 + ν) * (1 - 2 * ν) * Rᵢ^2 * -pressure(t) / (E * (Rₑ^2 - Rᵢ^2))
Avalue = A(t, Rᵢ, Rₑ, E, ν)
B(t, Rᵢ, Rₑ, E, ν) = (1 + ν) * Rᵢ^2 * Rₑ^2 * -pressure(t) / (E * (Rₑ^2 - Rᵢ^2))
Bvalue = B(t, Rᵢ, Rₑ, E, ν)
@show uᵣᵢ_analytic = Avalue * Rᵢ + Bvalue / Rᵢ
@show uᵣₑ_analytic = Avalue * Rₑ + Bvalue / Rₑ
nothing
