# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example  
#----------------------------------------------------
using LinearAlgebra: norm
using ONSAS
using Test: @test, @testset

## scalar parameters (dimensions in mm an MPa)
L = 75; # cylinder length in 𝐞ₖ mm
Rᵢ = 100; # inner radius in mm
Rₑ = 200; # outer radius in mm
p = 1e-1; # internal pressure in MPa
E = 210;  # Young modulus in MPa
ν = 0.3;  # Poisson ratio

# Run gmsh to generate the mesh
command = `gmsh -3 examples/cylinder_internal_pressure/cylinder.geo`;
run(command);

"Builds the `Structure`."
function cylinder_structure(E::Real, ν::Real, p::Real, L::Real, Rᵢ::Real, Rₑ::Real)
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
    bc₄ = LocalPressureBoundaryCondition([:u], t -> [p * t], bc₄_label)
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
structure = cylinder_structure(E, ν, p, L, Rᵢ, Rₑ);
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
states_lin_sol = solve!(linear_analysis);
# Select a face to test the boundary conditions 
"Tests boundary conditions and symmetry for a slice at x = Lᵢ."
function test_slice(Lᵢ::Real; atol::Real=ATOL)
    # Generic surface s at x = Lᵢ 
    rand_R = rand() * (Rₑ - Rᵢ) + Rᵢ
    rand_θ = rand() * 2 * π
    # Random point ∈ axis x
    prandᵢ = [rand_R, 0.0, Lᵢ]
    # Random point ∈ axis y
    prandⱼ = [0.0, rand_R, Lᵢ]
    # Random point between the internal and external surface
    prand = [rand_R * cos(rand_θ), rand_R * sin(rand_θ), Lᵢ]
    # Vector of points to test
    vec_points = [prandᵢ, prandⱼ, prand]
    #
    point_evaluator = PointEvalHandler(mesh(structure), vec_points)
    displacements_s = displacements(states_lin_sol, point_evaluator)
    # Check uₖ = 0 ∀ p ∈ s
    @test all([≈(norm(displacements_p[3]), 0.0, atol=atol) for displacements_p in displacements_s])
    # Check uᵢ = 0 ∀ p ∈ s & ∈ axis y
    displacements_axis_y_uᵢ = displacements(states_lin_sol, point_evaluator, 1)[1]
    @test all([≈(norm(displacements_p_axis_y_uᵢ), 0.0, atol=atol) for displacements_p_axis_y_uᵢ in displacements_axis_y_uᵢ])
    # Check uⱼ = 0 ∀ p ∈ s & ∈ axis x
    displacements_axis_y_uⱼ = displacements(states_lin_sol, point_evaluator, 2)[2]
    @test all([≈(norm(displacements_p_axis_y_uⱼ), 0.0, atol=atol) for displacements_p_axis_y_uⱼ in displacements_axis_y_uⱼ])
end
test_slice(L / 2, atol=ATOL)

