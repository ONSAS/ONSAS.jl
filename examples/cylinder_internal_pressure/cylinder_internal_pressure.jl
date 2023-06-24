# --------------------------------------------------
# Cylinder submitted to an Internal Pressure Example
#----------------------------------------------------
using LinearAlgebra, Test, Suppressor
using ONSAS

# Mesh with Gmsh.jl
include("cylinder_mesh.jl")

"Runs the cylinder with internal pressure example."
function run_cylinder_internal_pressure_example()
    ## scalar parameters (dimensions in mm an MPa)
    Lₖ = 30 # cylinder length in 𝐞ₖ mm
    Rᵢ = 100 # inner radius in mm
    Rₑ = 200 # outer radius in mm
    p = 10 # internal pressure in MPa
    E = 210.0  # Young modulus in MPa
    ν = 0.3  # Poisson ratio
    pressure(t::Real) = p * t
    ## number of steps
    NSTEPS = 9
    ## tolerances for testing
    ATOL = 1e-2 * (Rₑ - Rᵢ)
    ## Plot results
    PLOT_RESULTS = false
    ## Refinement mesh factor
    ms = 0.8 # ms = 2.5 for a Heavy one
    # -------------------------------
    # Structure
    # -------------------------------
    "Builds the `Structure`."
    function cylinder_structure(material::AbstractMaterial,
                                Lₖ::Real, Rᵢ::Real, Rₑ::Real,
                                pressure::Function; ms::Real)
        # -------------------------------
        # Physical entities labels
        # -------------------------------
        # material
        mat_label = label(material)
        # entities
        node_label = "node"
        faces_label = "triangle"
        elements_label = "tetrahedron"
        entities_labels = [node_label, faces_label, elements_label]
        # boundary conditions
        bc₁_label = "fixed-ui"
        bc₂_label = "fixed-uj"
        bc₃_label = "fixed-uk"
        bc₄_label = "pressure"
        bc_labels = [bc₁_label, bc₂_label, bc₃_label, bc₄_label]
        # mesh labels
        labels = [mat_label, entities_labels, bc_labels]
        # -------------------------------
        # Entities
        # -------------------------------
        # Entities types without assigned nodes, faces and elements
        vfaces = [TriangularFace(faces_label)]
        velems = [Tetrahedron(elements_label)]
        entities = StructuralEntity(velems, vfaces)
        # -------------------------------
        # Mesh
        # -------------------------------
        filename = "cylinder"
        local msh_path
        out = @capture_out begin
            msh_path = create_cylinder_mesh(Rᵢ, Rₑ, Lₖ, labels, filename, ms)
        end
        gmsh_println(out)
        msh_mesh = MshFile(msh_path)
        mesh = Mesh(msh_mesh, entities)
        # Dofs
        #--------------------------------
        dof_dim = 3
        dof_u_symbol = :u
        set_dofs!(mesh, dof_u_symbol, dof_dim)
        # -------------------------------
        # Boundary conditions
        # -------------------------------
        # Dirichlet boundary conditions
        bc₁ = FixedDof(:u, [1], bc₁_label)
        bc₂ = FixedDof(:u, [2], bc₂_label)
        bc₃ = FixedDof(:u, [3], bc₃_label)
        # Neumann boundary conditions
        bc₄ = Pressure(:u, pressure, bc₄_label)
        boundary_conditions = StructuralBoundaryCondition(bc₁, bc₂, bc₃, bc₄)
        # Assign boundary conditions to the ones defined in the mesh
        apply!(boundary_conditions, mesh)
        # -------------------------------
        # Materials
        # -------------------------------
        materials = StructuralMaterial(material)
        apply!(materials, mesh)
        # -------------------------------
        # Structure
        # -------------------------------
        Structure(mesh, materials, boundary_conditions)
    end
    # -------------------------------
    # Materials
    # -------------------------------
    mat_label = "mat"
    linear_material = IsotropicLinearElastic(E, ν, mat_label)
    cylinder = cylinder_structure(linear_material, Lₖ, Rᵢ, Rₑ, pressure; ms=ms)
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    "Defines an structural `AbstractStaticAnalysis`."
    function static_analysis(structure::Structure,
                             analysis::Type{<:AbstractStaticAnalysis};
                             NSTEPS::Int)
        analysis(structure; NSTEPS=NSTEPS)
    end
    # -----------------------------------------------
    # Case 1 - Static linear elastic case
    #------------------------------------------------
    linear_analysis = static_analysis(cylinder, LinearStaticAnalysis; NSTEPS=NSTEPS)
    # -------------------------------
    # Numerical solution
    # -------------------------------
    states_lin_sol = solve!(linear_analysis)
    # get time vector or load factors
    λᵥ = load_factors(linear_analysis)
    # Get the solution at a random point
    "Return a rand point in the cylinder (R, θ, L)."
    function rand_point_cylinder(Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ, Lₖ::Real=Lₖ)
        [rand() * (Rₑ - Rᵢ) + Rᵢ, rand() * 2 * π, rand() * Lₖ]
    end
    # Get the internal radial displacement at p = (0, Rᵢ, 0)
    cylinder_mesh = mesh(cylinder)
    nᵢ = nodes(cylinder_mesh)[4]
    uᵣ_numeric_nᵢ = displacements(states_lin_sol, nᵢ, 2)
    # Get the external radial displacement at p = (-Rₑ, 0, Lₖ)
    nₑ = nodes(cylinder_mesh)[15]
    uᵣ_numeric_nₑ = displacements(states_lin_sol, nₑ, 1)
    # Generate a random point
    rand_R, rand_θ, rand_z = rand_point_cylinder()
    # rand_R = Rᵢ;
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
                E::Real=E, ν::Real=ν, p::Function=pressure,
                Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ)
        "Constant A for the analytic solution."
        function A(t::Real, Rᵢ::Real, Rₑ::Real, E::Real, ν::Real, p::Function)
            (1 + ν) * (1 - 2 * ν) * Rᵢ^2 * p(t) / (E * (Rₑ^2 - Rᵢ^2))
        end
        "Constant B for the analytic solution."
        function B(t::Real, Rᵢ::Real, Rₑ::Real, E::Real, ν::Real, p::Function)
            (1 + ν) * Rᵢ^2 * Rₑ^2 * p(t) / (E * (Rₑ^2 - Rᵢ^2))
        end
        A(t, Rᵢ, Rₑ, E, ν, p) * r + B(t, Rᵢ, Rₑ, E, ν, p) / r
    end
    uᵣ_analytic_nᵢ = [uᵣ(Rᵢ, λᵢ) for λᵢ in λᵥ]
    uᵣ_analytic_nₑ = [uᵣ(Rₑ, λᵢ) for λᵢ in λᵥ]
    uᵣ_analytic_p_rand = [uᵣ(rand_R, λᵢ) for λᵢ in λᵥ]
    #-----------------------------
    # Test booleans - Case 1
    #-----------------------------
    # Test symmetry and boundary conditions for a random slice
    #-------------------------------------------
    function test_solution_at_slice(sol::AbstractSolution=states_lin_sol;
                                    atol::Real=ATOL, atolr=ATOLR,
                                    Rᵢ::Real=Rᵢ, Rₑ::Real=Rₑ, Lₖ::Real=Lₖ)
        structure = ONSAS.structure(analysis(sol))
        # Generic surface s at z = Lₖ
        rand_R, rand_θ₁, Lₖ = rand_point_cylinder(Rᵢ, Rₑ, Lₖ)
        # Set by force Lₖ
        rand_θ₂ = rand() * 2 * π
        # Random point ∈ axis x
        p_rand_in_axis_x = [rand_R, 0.0, Lₖ]
        # Random point ∈ axis y
        p_rand_in_axis_y = [0.0, rand_R, Lₖ]
        # Random point between the internal and external surface
        p_rand₁ = [rand_R * cos(rand_θ₁), rand_R * sin(rand_θ₁), Lₖ]
        p_rand₂ = [rand_R * cos(rand_θ₂), rand_R * sin(rand_θ₂), Lₖ]
        # Vector of points to test
        vec_points = [p_rand_in_axis_x, p_rand_in_axis_y, p_rand₁, p_rand₂]
        #
        point_evaluator = PointEvalHandler(mesh(structure), vec_points)
        U = displacements(sol, point_evaluator)
        # Check uₖ = 0 ∀ p ∈ s
        zero_uₖ = all([≈(norm(u[3]), 0.0; atol=atol) for u in U])
        # Check uᵢ = 0 ∀ p ∈ s & ∈ axis y
        index_p_rand_in_axis_y = findall([p == p_rand_in_axis_y for p in vec_points])
        Uᵢ_in_axis_y = getindex(displacements(sol, point_evaluator, 1), index_p_rand_in_axis_y)
        zero_uᵢ_axis_y = all([≈(norm(uᵢ_p_in_axis_y), 0.0; atol=atol)
                              for uᵢ_p_in_axis_y in Uᵢ_in_axis_y])
        # Check uⱼ = 0 ∀ p ∈ s & ∈ axis x
        index_p_rand_in_axis_x = findall([p == p_rand_in_axis_x for p in vec_points])
        Uⱼ_in_axis_x = getindex(displacements(sol, point_evaluator, 2), index_p_rand_in_axis_x)
        zero_uⱼ_axis_x = all([≈(norm(uⱼ_p_in_axis_y), 0.0; atol=atol)
                              for uⱼ_p_in_axis_y in Uⱼ_in_axis_x])
        # Check uᵣ(r,θ₁) =  uᵣ(r,θ₁)  at last time
        rand₁_index = 3
        uᵣ_rand₁ = sum(last.(U[rand₁_index][1:2]) .^ 2)
        rand₂_index = 4
        uᵣ_rand₂ = sum(last.(U[rand₂_index][1:2]) .^ 2)
        uᵣ_not_depends_on_θ = ≈(uᵣ_rand₁, uᵣ_rand₂; atol=atolr)

        uᵣ_not_depends_on_θ, zero_uₖ, zero_uᵢ_axis_y, zero_uⱼ_axis_x
    end
    # Test symmetry and boundary conditions
    test_bools_symmetry_linear = test_solution_at_slice(states_lin_sol; atol=ATOL, atolr=10 * ATOL)
    uᵣ_not_depends_on_θ_linear, zero_uₖ_linear, zero_uᵢ_axis_y_linear, zero_uⱼ_axis_x_linear = test_bools_symmetry_linear
    # -----------------------------------------------
    # Case 2 - Static non-linear elastic case
    #----------------------------------------------
    svk_material = Svk(; E=E, ν=ν, label=mat_label)
    replace!(cylinder, svk_material)
    nonlinear_analysis = static_analysis(cylinder, NonLinearStaticAnalysis; NSTEPS=NSTEPS)
    #  Non-linear analysis
    # -------------------------------
    tols = ConvergenceSettings(; rel_U_tol=1e-8, rel_res_force_tol=1e-8, max_iter=30)
    alg = NewtonRaphson(tols)
    states_nonlinear_sol = solve!(nonlinear_analysis, alg)
    # Get the internal radial displacement at p = (0, Rᵢ, 0)
    uᵣ_numeric_nonlinear_nᵢ = displacements(states_nonlinear_sol, nᵢ, 2)
    # Get the external radial displacement at p = (-Rₑ, 0, Lₖ)
    uᵣ_numeric_nonlinear_nₑ = displacements(states_nonlinear_sol, nₑ, 1)
    # Test symmetry and boundary conditions
    test_silce_bools = test_solution_at_slice(states_nonlinear_sol; atol=ATOL, atolr=10 * ATOL)
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
    fig = plot(vec_p, uᵣ_numeric_nᵢ; label="numeric linear uᵣ n=(0, Rᵢ, 0)",
               legend=:topleft, color=:orange, lw=2, ls=:dash, markershape=:circle)
    plot!(fig,
          vec_p, -uᵣ_numeric_nₑ; label="numeric linear uᵣ n=(-Rₑ, 0 , Lₖ)",
          legend=:topleft, color=:skyblue, lw=2, ls=:solid, markershape=:square)
    plot!(fig,
          vec_p, uᵣ_analytic_nᵢ; label="analytic linear uᵣ(Rᵢ)",
          legend=:topleft, color=:black, lw=2, ls=:dash, markershape=:none)
    plot!(fig,
          vec_p, uᵣ_analytic_nₑ; label="analytic linear uᵣ(Rₑ)",
          legend=:topleft, color=:black, lw=2, ls=:solid)
    # Plot comparing linear and non linear solutions
    plot!(fig,
          vec_p_non_in, uᵣ_numeric_nonlinear_nᵢ; label="non-linear uᵣ(0, Rᵢ, 0)",
          color=:red, lw=2, marker=:circle, markersize=3)
    plot!(fig,
          vec_p_non_in, -uᵣ_numeric_nonlinear_nₑ; label="non-linear uᵣ(-Rₑ, 0 , Lₖ)",
          color=:blue, lw=2, marker=:circle, markersize=3)
    # add labels
    xlabel!("λᵥ [MPa]")
    ylabel!("uᵣ [mm]")
    display(fig)
end

run_cylinder_internal_pressure_example()
