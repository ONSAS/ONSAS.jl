# --------------------------------------------------------------------------
# Uniaxial Extension ExampleExercise 4 from section 6.5 in (Holzapfel,2000).
# For notation see: https://onsas.github.io/ONSAS.m/dev/examples/uniaxialExtension/
# --------------------------------------------------------------------------
using Test, LinearAlgebra, Suppressor, Roots
using ONSAS

include("uniaxial_mesh.jl") # Mesh Cube with Gmsh.jl

function run_uniaxial_extension()
	## scalar parameters
	E = 1.0                    # Young modulus in Pa
	ν = 0.3                    # Poisson's ratio
	p = 3                      # Tension load in Pa
	Lᵢ = 2.0                   # Dimension in x of the box in m
	Lⱼ = 1.0                   # Dimension in y of the box in m
	Lₖ = 1.0                   # Dimension in z of the box in m
	ms = 0.5                   # Refinement factor for the mesh
	RTOL = 1e-4                # Relative tolerance for tests
	# -----------------------------------------------
	# Case 1 - Manufactured mesh and `SVK` material
	#------------------------------------------------
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
	s₁_mesh = Mesh(; nodes = vec_nodes)
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
	append!(faces(s₁_mesh), vec_faces)
	## Entities
	t₁ = Tetrahedron(n₁, n₄, n₂, n₆, "tetra_1")
	t₂ = Tetrahedron(n₆, n₂, n₃, n₄, "tetra_2")
	t₃ = Tetrahedron(n₄, n₃, n₆, n₇, "tetra_3")
	t₄ = Tetrahedron(n₄, n₁, n₅, n₆, "tetra_4")
	t₅ = Tetrahedron(n₄, n₆, n₅, n₈, "tetra_5")
	t₆ = Tetrahedron(n₄, n₇, n₆, n₈, "tetra_6")
	vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
	append!(elements(s₁_mesh), vec_elems)
	# -------------------------------
	# Dofs
	#--------------------------------
	dof_dim = 3
	dof_u_symbol = :u
	set_dofs!(s₁_mesh, dof_u_symbol, dof_dim)
	# -------------------------------
	# Materials
	# -------------------------------
	svk = SVK(; E = E, ν = ν, label = "svk")
	s₁_materials = StructuralMaterial(svk => [t₁, t₂, t₃, t₄, t₅, t₆])
	# -------------------------------
	# Boundary conditions
	# -------------------------------
	# Fixed dofs
	bc₁_label = "fixed-ux"
	bc₁ = FixedDof(:u, [1], bc₁_label)
	bc₂_label = "fixed-uj"
	bc₂ = FixedDof(:u, [2], bc₂_label)
	bc₃_label = "fixed-uk"
	bc₃ = FixedDof(:u, [3], bc₃_label)
	# Load
	bc₄_label = "tension"
	bc₄ = GlobalLoad(:u, t -> [p * t, 0, 0], bc₄_label)
	# Crete boundary conditions struct
	s₁_boundary_conditions = StructuralBoundaryCondition(bc₁ => [f₃, f₄], bc₂ => [f₅, f₆],
		bc₃ => [f₇, f₈], bc₄ => [f₁, f₂])
	bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
	# -------------------------------
	# Structure
	# -------------------------------
	s₁ = Structure(s₁_mesh, s₁_materials, s₁_boundary_conditions)
	# -------------------------------
	# Structural Analysis
	# -------------------------------
	# Final load factor
	NSTEPS = 8
	sa₁ = NonLinearStaticAnalysis(s₁; NSTEPS = NSTEPS)
	# -------------------------------
	# Algorithm
	# -------------------------------
	tol_f = 1e-8
	tol_u = 1e-8
	max_iter = 30
	tols = ConvergenceSettings(tol_u, tol_f, max_iter)
	nr = NewtonRaphson(tols)
	# -------------------------------
	# Numerical solution
	# -------------------------------
	states_sol_case₁ = ONSAS.solve!(sa₁, nr)

	# -------------------------------
	# Write vtk files
	# -------------------------------
	vtk_path = joinpath("examples", "uniaxial_extension", "uniaxial")
	write_vtk(states_sol_case₁, vtk_path)
	"Computes numeric solution α, β and γ for analytic validation."
	function αβγ_numeric(states_sol::AbstractSolution)
		s = structure(analysis(states_sol))
		# Node at (Lᵢ, Lⱼ, Lₖ)
		n₇ = nodes(s)[7]
		numerical_uᵢ = displacements(states_sol_case₁, n₇, 1)
		numerical_α = 1 .+ numerical_uᵢ / Lᵢ
		# Displacements in the y (component 2) axis at node 7
		numerical_uⱼ = displacements(states_sol_case₁, n₇, 2)
		numerical_β = 1 .+ numerical_uⱼ / Lⱼ
		# Displacements in the z (component 3) axis at node 7
		numerical_uₖ = displacements(states_sol_case₁, n₇, 3)
		numerical_γ = 1 .+ numerical_uₖ / Lₖ
		numerical_α, numerical_β, numerical_γ, numerical_uᵢ, numerical_uⱼ, numerical_uₖ
	end
	# Numeric solution for testing
	numeric_α_case₁, numeric_β_case₁, numeric_γ_case₁, numeric_uᵢ_case₁, _, _ = αβγ_numeric(states_sol_case₁)
	# Extract ℙ and ℂ from the last state using a random element
	e = rand(elements(s₁))
	# Cosserat or second Piola-Kirchhoff stress tensor
	ℙ_numeric_case₁ = last(stress(states_sol_case₁, e))
	# Right hand Cauchy strain tensor
	ℂ_numeric_case₁ = last(strain(states_sol_case₁, e))
	# Load factors
	numeric_λᵥ_case₁ = load_factors(sa₁)

	# -----------------------------------------------
	# Case 2 - GMSH mesh and `HyperElastic` material
	#------------------------------------------------
	# -------------------------------
	# Materials
	# -------------------------------
	# Define a new HyperElastic material from the strain energy function
	strain_energy_svk(𝔼::AbstractMatrix, λ::Real, G::Real) = (λ / 2) * tr(𝔼)^2 + G * tr(𝔼^2)
	λ, G = lame_parameters(svk)
	params = [λ, G] # The order must be the same defined in the strain energy(splatting)
	mat_label = "svkHyper"
	svk_hyper_elastic = HyperElastic(params, strain_energy_svk, mat_label)
	# Material types without assigned elements
	s_materials = StructuralMaterial(svk_hyper_elastic)
	# -------------------------------
	# Boundary Conditions
	# -------------------------------
	# Redefine the load boundary condition
	bc₄ = Pressure(:u, t -> -p * t, bc₄_label)
	# BoundaryConditions types without assigned node, feces and elements
	s_boundary_conditions = StructuralBoundaryCondition(bc₁, bc₂, bc₃, bc₄)
	# -------------------------------
	# Entities
	# -------------------------------
	# Entities types without assigned nodes, faces and elements
	faces_label = "triangle"
	elems_label = "tetrahedron"
	vfaces = [TriangularFace(faces_label)]
	velems = [Tetrahedron(elems_label)]
	s_entities = StructuralEntity(velems, vfaces)
	entities_labels = [faces_label, elems_label]
	# -------------------------------
	# Mesh
	# -------------------------------
	filename = "uniaxial_extension"
	labels = [mat_label, entities_labels, bc_labels]
	local mesh_path
	output = @capture_out begin
		mesh_path = create_uniaxial_mesh(Lᵢ, Lⱼ, Lₖ, labels, filename, ms)
	end
	gmsh_println(output)
	msh_file = MshFile(mesh_path)
	# -------------------------------
	# Structure
	# -------------------------------
	s₂ = Structure(msh_file, s_materials, s_boundary_conditions, s_entities)
	# Final load factor
	sa₂ = NonLinearStaticAnalysis(s₂; NSTEPS = NSTEPS)
	# -------------------------------
	# Numerical solution
	# -------------------------------
	states_sol_case₂ = ONSAS.solve!(sa₂, nr)
	# Numeric solution for testing
	numeric_α_case₂, numeric_β_case₂, numeric_γ_case₂, numeric_uᵢ_case₂, _, _ = αβγ_numeric(states_sol_case₂)
	# Extract ℙ and ℂ from the last state using a random element
	e = rand(elements(s₂))
	# Cosserat or second Piola-Kirchhoff stress tensor
	ℙ_numeric_case₂ = last(stress(states_sol_case₂, e))
	# Right hand Cauchy strain tensor
	ℂ_numeric_case₂ = last(strain(states_sol_case₂, e))
	# Load factors
	numeric_λᵥ_case₂ = load_factors(sa₂)
	#-----------------------------
	# Analytic solution
	#-----------------------------
	"Computes displacements numeric solution uᵢ, uⱼ and uₖ for analytic validation."
	function u_ijk_numeric(numerical_α::Vector{<:Real}, numerical_β::Vector{<:Real},
		numerical_γ::Vector{<:Real},
		x::Real, y::Real, z::Real)
		x * (numerical_α .- 1), y * (numerical_β .- 1), z * (numerical_γ .- 1)
	end
	# Test with load factors
	"Analytic load factor solution for the displacement `uᵢ` towards `x` axis at node `n₆`."
	function load_factors_analytic(uᵢ::Real, p::Real = p, E::Real = E, Lᵢ::Real = Lᵢ)
		1 / p * E * 0.5 * ((1 + uᵢ / Lᵢ)^3 - (1 + uᵢ / Lᵢ))
	end
	# Compute load factors with numerical solutions
	analytics_λᵥ_case₁ = load_factors_analytic.(numeric_uᵢ_case₁)
	analytics_λᵥ_case₂ = load_factors_analytic.(numeric_uᵢ_case₂)
	# Test last step σ and ϵ
	@test load_factors(sa₁) == load_factors(sa₂)
	α_analytic = find_zero(α -> E / 2 * α * (α^2 - 1) - p * last(load_factors(sa₁)), 1e-2)
	β_analytic = sqrt(-ν * (α_analytic^2 - 1) + 1)
	# Gradient tensor
	# 𝑢 = (αx, βy, γz)
	𝔽_analytic = [α_analytic 0 0; 0 β_analytic 0; 0 0 β_analytic]
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
	# -------------------------------
	# Interpolator tests for Case 2
	#--------------------------------
	rand_point = [[rand() * Lᵢ, rand() * Lⱼ, rand() * Lₖ]]
	eval_handler_rand = PointEvalHandler(mesh(s₂), rand_point)
	# Compute analytic solution at a random point
	uᵢ_case₂, uⱼ_case₂, uₖ_case₂ = u_ijk_numeric(numeric_α_case₂, numeric_β_case₂, numeric_γ_case₂,
		rand_point[]...)
	rand_point_uᵢ = displacements(states_sol_case₂, eval_handler_rand, 1)
	rand_point_uⱼ = displacements(states_sol_case₂, eval_handler_rand, 2)
	rand_point_uₖ = displacements(states_sol_case₂, eval_handler_rand, 3)
	#-----------------------------
	# Test booleans
	#-----------------------------
	@testset "Case 1 Uniaxial Extension Example" begin
		@test analytics_λᵥ_case₁ ≈ numeric_λᵥ_case₁ rtol = RTOL
		@test ℙ_analytic ≈ ℙ_numeric_case₁ rtol = RTOL
		@test α_analytic ≈ last(numeric_α_case₁) rtol = RTOL
		@test ℂ_analytic ≈ ℂ_numeric_case₁ rtol = RTOL
	end

	@testset "Case 2 Uniaxial Extension Example" begin
		@test analytics_λᵥ_case₂ ≈ numeric_λᵥ_case₂ rtol = RTOL
		@test ℂ_analytic ≈ ℂ_numeric_case₂ rtol = RTOL
		@test ℙ_analytic ≈ ℙ_numeric_case₂ rtol = RTOL
		@test β_analytic ≈ last(numeric_β_case₂) rtol = RTOL
		# Interpolations
		@test uᵢ_case₂ ≈ rand_point_uᵢ rtol = RTOL
		@test uⱼ_case₂ ≈ rand_point_uⱼ rtol = RTOL
		@test uₖ_case₂ ≈ rand_point_uₖ rtol = RTOL
	end
end

run_uniaxial_extension()
