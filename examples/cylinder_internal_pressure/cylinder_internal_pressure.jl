# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example
#----------------------------------------------------
using LinearAlgebra, Test, Suppressor
using ONSAS

# Mesh with Gmsh.jl
include("cylinder_mesh.jl")

"Return problem parameters"
function parameters()
	## scalar parameters (dimensions in mm an MPa)
	Lz = 30                 # cylinder length in z mm
	Ri = 100                # inner radius in mm
	Re = 200                # outer radius in mm
	p = 10                  # internal pressure in MPa
	pressure(t::Real) = p * t
	material_label = "mat"  # Material label used to replace the material structure
	E = 210                 # Young modulus in MPa
	ν = 0.3                 # Poisson ratio
	NSTEPS = 9              # Number of load steps
	ATOL = 1e-2 * (Re - Ri) # Absolut tolerances for testing
	ms = 0.8                # ms = 2.5 for much more refined mesh (approx 200.000 elems)
	PLOT_RESULTS = false    # boolean to plot results
	(; Lz, Ri, Re, p, ν, E, pressure, ATOL, NSTEPS, PLOT_RESULTS, ms, material_label)
end;

#= -----------------------------------------------------------
Two cases are considered:
Case 1 - Linear Elastic Analysis with Analytic Soltuion
Case 2 - Hyper Elastic Analysis to check the solution's symmetry
-------------------------------------------------------------=#
abstract type AbstractCase end
struct FirstCase <: AbstractCase end
struct SecondCase <: AbstractCase end


"Return the problem structural model"
function structure(material::AbstractMaterial)
	(; Lz, Ri, Re, pressure, ms) = parameters()
	# -------------------------------
	# Materials
	# -------------------------------
	mat_label = label(material)
	materials = StructuralMaterial(material)
	# -------------------------------
	# Boundary conditions
	# -------------------------------
	# Dirichlet boundary conditions
	bc1_label = "fixed-ui"
	bc1 = FixedDof(:u, [1], bc1_label)
	bc2_label = "fixed-uj"
	bc2 = FixedDof(:u, [2], bc2_label)
	bc3_label = "fixed-uk"
	bc3 = FixedDof(:u, [3], bc3_label)
	# Neumann boundary conditions
	bc4_label = "pressure"
	bc4 = Pressure(:u, pressure, bc4_label)
	boundary_conditions = StructuralBoundaryCondition(bc1, bc2, bc3, bc4)
	bc_labels = [bc1_label, bc2_label, bc3_label, bc4_label]
	# -------------------------------
	# Entities
	# -------------------------------
	# Entities types without assigned nodes, faces and elements
	node_label = "node"
	faces_label = "triangle"
	elements_label = "tetrahedron"
	vfaces = [TriangularFace(faces_label)]
	velems = [Tetrahedron(elements_label)]
	entities_labels = [node_label, faces_label, elements_label]
	entities = StructuralEntity(velems, vfaces)
	# -------------------------------
	# Mesh
	# -------------------------------
	labels = [mat_label, entities_labels, bc_labels]
	filename = "cylinder"
	local msh_path
	out = @capture_out begin
		msh_path = create_cylinder_mesh(Ri, Re, Lz, labels, filename, ms)
	end
	gmsh_println(out)
	msh_mesh = MshFile(msh_path)
	mesh = Mesh(msh_mesh, entities)
	#--------------------------------
	# Dofs
	#--------------------------------
	dof_dim = 3
	dof_u_symbol = :u
	set_dofs!(mesh, dof_u_symbol, dof_dim)
	# -------------------------------
	# Structure
	# -------------------------------
	apply!(materials, mesh)
	apply!(boundary_conditions, mesh)

	Structure(mesh, materials, boundary_conditions)
end;


"Return the problem solution"
function solve(::FirstCase)
	(; NSTEPS, E, ν, material_label) = parameters()
	liner_elastic = IsotropicLinearElastic(E, ν, material_label)
	# -------------------------------
	# Structural Analysis
	# -------------------------------
	s = structure(liner_elastic)
	sa = LinearStaticAnalysis(s; NSTEPS)
	# -------------------------------
	# Numerical solution
	# -------------------------------
	solve!(sa)
end;


"Return the problem solution"
function solve(::SecondCase)
	(; NSTEPS, E, ν, material_label) = parameters()
	# The material is replaced just to test the replace method
	liner_elastic = IsotropicLinearElastic(E, ν, material_label)
	s = structure(liner_elastic)
	svk_material = SVK(; E = E, ν = ν, label = material_label)
	replace!(s, svk_material)
	# -------------------------------
	# Structural Analysis
	# -------------------------------
	sa = NonLinearStaticAnalysis(s; NSTEPS)
	# -------------------------------
	# Solver
	# -------------------------------
	tol_f = 1e-10
	tol_u = 1e-10
	max_iter = 30
	tols = ConvergenceSettings(tol_u, tol_f, max_iter)
	nr = NewtonRaphson(tols)
	# -------------------------------
	# Numerical solution
	# -------------------------------
	solve!(sa, nr)
end;





"Runs the cylinder with internal pressure example."
function run_cylinder_internal_pressure_example()
	# -------------------------------
	# Structure
	# -------------------------------

	# -------------------------------
	# Materials
	# -------------------------------
	linear_material = IsotropicLinearElastic(E, ν, mat_label)
	cylinder = cylinder_structure(linear_material, Lz, Ri, Re, pressure; ms = ms)
	# -------------------------------
	# Structural Analysis
	# -------------------------------
	"Defines an structural `AbstractStaticAnalysis`."
	function static_analysis(structure::Structure,
		analysis::Type{<:AbstractStaticAnalysis};
		NSTEPS::Int)
		analysis(structure; NSTEPS = NSTEPS)
	end
	# -----------------------------------------------
	# Case 1 - Static linear elastic case
	#------------------------------------------------
	linear_analysis = static_analysis(cylinder, LinearStaticAnalysis; NSTEPS = NSTEPS)
	# -------------------------------
	# Numerical solution
	# -------------------------------
	states_lin_sol = solve!(linear_analysis)
	# get time vector or load factors
	λᵥ = load_factors(linear_analysis)
	# Get the solution at a random point
	"Return a rand point in the cylinder (R, θ, L)."
	function rand_point_cylinder(Ri::Real = Ri, Re::Real = Re, Lz::Real = Lz)
		[rand() * (Re - Ri) + Ri, rand() * 2 * π, rand() * Lz]
	end
	# Get the internal radial displacement at p = (0, Ri, 0)
	cylinder_mesh = mesh(cylinder)
	nᵢ = nodes(cylinder_mesh)[4]
	uᵣ_numeric_nᵢ = displacements(states_lin_sol, nᵢ, 2)
	# Get the external radial displacement at p = (-Re, 0, Lz)
	nₑ = nodes(cylinder_mesh)[15]
	uᵣ_numeric_nₑ = displacements(states_lin_sol, nₑ, 1)
	# Generate a random point
	rand_R, rand_θ, rand_z = rand_point_cylinder()
	# rand_R = Ri;
	# rand_θ = 0.0
	p_rand = Point(rand_R * cos(rand_θ), rand_R * sin(rand_θ), rand_z)
	# Displacements at p
	point_evaluator = PointEvalHandler(cylinder_mesh, p_rand)
	uᵢ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 1)
	uⱼ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 2)
	uₖ_numeric_p_rand = displacements(states_lin_sol, point_evaluator, 3)
	#
	uᵣ_numeric_p_rand = sqrt.(@. uᵢ_numeric_p_rand^2 + uⱼ_numeric_p_rand^2)
	# -------------------------------
	# Analytic solution
	# -------------------------------
	"Analytic radial displacements uᵣ at radius`r` and time `t`."
	function uᵣ(r::Real, t::Real,
		E::Real = E, ν::Real = ν, p::Function = pressure,
		Ri::Real = Ri, Re::Real = Re)
		"Constant A for the analytic solution."
		function A(t::Real, Ri::Real, Re::Real, E::Real, ν::Real, p::Function)
			(1 + ν) * (1 - 2 * ν) * Ri^2 * p(t) / (E * (Re^2 - Ri^2))
		end
		"Constant B for the analytic solution."
		function B(t::Real, Ri::Real, Re::Real, E::Real, ν::Real, p::Function)
			(1 + ν) * Ri^2 * Re^2 * p(t) / (E * (Re^2 - Ri^2))
		end
		A(t, Ri, Re, E, ν, p) * r + B(t, Ri, Re, E, ν, p) / r
	end
	uᵣ_analytic_nᵢ = [uᵣ(Ri, λᵢ) for λᵢ in λᵥ]
	uᵣ_analytic_nₑ = [uᵣ(Re, λᵢ) for λᵢ in λᵥ]
	uᵣ_analytic_p_rand = [uᵣ(rand_R, λᵢ) for λᵢ in λᵥ]
	#-----------------------------
	# Test booleans - Case 1
	#-----------------------------
	# Test symmetry and boundary conditions for a random slice
	#-------------------------------------------
	function test_solution_at_slice(sol::AbstractSolution = states_lin_sol;
		atol::Real = ATOL, atolr = ATOLR,
		Ri::Real = Ri, Re::Real = Re, Lz::Real = Lz)
		structure = ONSAS.structure(analysis(sol))
		# Generic surface s at z = Lz
		rand_R, rand_θ₁, Lz = rand_point_cylinder(Ri, Re, Lz)
		# Set by force Lz
		rand_θ₂ = rand() * 2 * π
		# Random point ∈ axis x
		p_rand_in_axis_x = [rand_R, 0.0, Lz]
		# Random point ∈ axis y
		p_rand_in_axis_y = [0.0, rand_R, Lz]
		# Random point between the internal and external surface
		p_rand₁ = [rand_R * cos(rand_θ₁), rand_R * sin(rand_θ₁), Lz]
		p_rand₂ = [rand_R * cos(rand_θ₂), rand_R * sin(rand_θ₂), Lz]
		# Vector of points to test
		vec_points = [p_rand_in_axis_x, p_rand_in_axis_y, p_rand₁, p_rand₂]
		#
		point_evaluator = PointEvalHandler(mesh(structure), vec_points)
		U = displacements(sol, point_evaluator)
		# Check uₖ = 0 ∀ p ∈ s
		zero_uₖ = all([≈(norm(u[3]), 0.0; atol = atol) for u in U])
		# Check uᵢ = 0 ∀ p ∈ s & ∈ axis y
		index_p_rand_in_axis_y = findall([p == p_rand_in_axis_y for p in vec_points])
		Uᵢ_in_axis_y = getindex(displacements(sol, point_evaluator, 1), index_p_rand_in_axis_y)
		zero_uᵢ_axis_y = all([≈(norm(uᵢ_p_in_axis_y), 0.0; atol = atol)
							  for uᵢ_p_in_axis_y in Uᵢ_in_axis_y])
		# Check uⱼ = 0 ∀ p ∈ s & ∈ axis x
		index_p_rand_in_axis_x = findall([p == p_rand_in_axis_x for p in vec_points])
		Uⱼ_in_axis_x = getindex(displacements(sol, point_evaluator, 2), index_p_rand_in_axis_x)
		zero_uⱼ_axis_x = all([≈(norm(uⱼ_p_in_axis_y), 0.0; atol = atol)
							  for uⱼ_p_in_axis_y in Uⱼ_in_axis_x])
		# Check uᵣ(r,θ₁) =  uᵣ(r,θ₁)  at last time
		rand₁_index = 3
		uᵣ_rand₁ = sum(last.(U[rand₁_index][1:2]) .^ 2)
		rand₂_index = 4
		uᵣ_rand₂ = sum(last.(U[rand₂_index][1:2]) .^ 2)
		uᵣ_not_depends_on_θ = ≈(uᵣ_rand₁, uᵣ_rand₂; atol = atolr)

		uᵣ_not_depends_on_θ, zero_uₖ, zero_uᵢ_axis_y, zero_uⱼ_axis_x
	end
	# Test symmetry and boundary conditions
	test_bools_symmetry_linear = test_solution_at_slice(states_lin_sol; atol = ATOL, atolr = 10 * ATOL)
	uᵣ_not_depends_on_θ_linear, zero_uₖ_linear, zero_uᵢ_axis_y_linear, zero_uⱼ_axis_x_linear = test_bools_symmetry_linear
	# -----------------------------------------------
	# Case 2 - Static non-linear elastic case
	#----------------------------------------------
	svk_material = SVK(; E = E, ν = ν, label = mat_label)
	replace!(cylinder, svk_material)
	nonlinear_analysis = static_analysis(cylinder, NonLinearStaticAnalysis; NSTEPS = NSTEPS)
	#  Non-linear analysis
	# -------------------------------
	tols = ConvergenceSettings(; rel_U_tol = 1e-8, rel_res_force_tol = 1e-8, max_iter = 30)
	alg = NewtonRaphson(tols)
	states_nonlinear_sol = solve!(nonlinear_analysis, alg)
	# Get the internal radial displacement at p = (0, Ri, 0)
	uᵣ_numeric_nonlinear_nᵢ = displacements(states_nonlinear_sol, nᵢ, 2)
	# Get the external radial displacement at p = (-Re, 0, Lz)
	uᵣ_numeric_nonlinear_nₑ = displacements(states_nonlinear_sol, nₑ, 1)
	# Test symmetry and boundary conditions
	test_silce_bools = test_solution_at_slice(states_nonlinear_sol; atol = ATOL, atolr = 10 * ATOL)
	uᵣ_not_depends_on_θ_nonlinear, zero_uₖ_case2_nonlinear, zero_uᵢ_axis_y_nonlinear, zero_uⱼ_axis_x_nonlinear = test_silce_bools
	#-----------------------------
	# Test booleans
	#-----------------------------
	@testset "Case 1: Linear Analysis " begin
		@test uᵣ_not_depends_on_θ_linear
		@test zero_uₖ_linear
		@test zero_uᵢ_axis_y_linear
		@test zero_uⱼ_axis_x_linear
		@test uᵣ_numeric_p_rand ≈ uᵣ_analytic_p_rand atol = ATOL
		@test uᵣ_analytic_nᵢ ≈ uᵣ_numeric_nᵢ atol = ATOL
		@test uᵣ_analytic_nₑ ≈ -uᵣ_numeric_nₑ atol = ATOL
	end
	@testset "Case 2: Non-Linear Analysis " begin
		@test uᵣ_not_depends_on_θ_nonlinear
		@test zero_uₖ_case2_nonlinear
		@test zero_uᵢ_axis_y_nonlinear
		@test zero_uⱼ_axis_x_nonlinear
	end
	#-----------------------------
	# Plot & plots
	#-----------------------------
	PLOT_RESULTS && plot_results(λᵥ, nonlinear_analysis, uᵣ_numeric_nᵢ, uᵣ_numeric_nₑ,
		uᵣ_analytic_nᵢ, uᵣ_analytic_nₑ)
end

"Plot cylinder with internal pressure results uᵣ vs λ"
function plot_results(λᵥ, nonlinear_analysis, uᵣ_numeric_nᵢ, uᵣ_numeric_nₑ, uᵣ_analytic_nᵢ,
	uᵣ_analytic_nₑ)
	vec_p = [pressure(λ) for λ in λᵥ]
	vec_p_non_in = [pressure(λ) for λ in load_factors(nonlinear_analysis)]
	fig = plot(vec_p, uᵣ_numeric_nᵢ; label = "numeric linear uᵣ n=(0, Ri, 0)",
		legend = :topleft, color = :orange, lw = 2, ls = :dash, markershape = :circle)
	plot!(fig,
		vec_p, -uᵣ_numeric_nₑ; label = "numeric linear uᵣ n=(-Re, 0 , Lz)",
		legend = :topleft, color = :skyblue, lw = 2, ls = :solid, markershape = :square)
	plot!(fig,
		vec_p, uᵣ_analytic_nᵢ; label = "analytic linear uᵣ(Ri)",
		legend = :topleft, color = :black, lw = 2, ls = :dash, markershape = :none)
	plot!(fig,
		vec_p, uᵣ_analytic_nₑ; label = "analytic linear uᵣ(Re)",
		legend = :topleft, color = :black, lw = 2, ls = :solid)
	# Plot comparing linear and non linear solutions
	plot!(fig,
		vec_p_non_in, uᵣ_numeric_nonlinear_nᵢ; label = "non-linear uᵣ(0, Ri, 0)",
		color = :red, lw = 2, marker = :circle, markersize = 3)
	plot!(fig,
		vec_p_non_in, -uᵣ_numeric_nonlinear_nₑ; label = "non-linear uᵣ(-Re, 0 , Lz)",
		color = :blue, lw = 2, marker = :circle, markersize = 3)
	# add labels
	xlabel!("λᵥ [MPa]")
	ylabel!("uᵣ [mm]")
	display(fig)
end

run_cylinder_internal_pressure_example()
