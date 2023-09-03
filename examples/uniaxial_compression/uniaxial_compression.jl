# ------------------------
# Von Misses Truss Example
# ------------------------
using Test, LinearAlgebra, Suppressor
using ONSAS

# Mesh Cube with Gmsh.jl
include(joinpath("..", "uniaxial_extension", "uniaxial_mesh.jl"))

"Return problem parameters"
function parameters()
    E = 1.0                     # Young modulus in Pa
    ν = 0.3                     # Poisson's ratio
    μ = G = E / (2 * (1 + ν))   # Second Lamé parameter
    K = E / (3 * (1 - 2 * ν))   # Bulk modulus
    p = 1                       # Pressure load in Pa
    Li = 2.0                    # Dimension in x of the box in m
    Lj = 1.0                    # Dimension in y of the box in m
    Lk = 1.0                    # Dimension in z of the box in m
    ms = 0.5                    # Refinement factor for the mesh
    RTOL = 1e-4                 # Relative tolerance for tests
    ATOL = 1e-10                # Absolute tolerance for tests
    NSTEPS = 9                  # Number of steps for the test
    (; μ, G, K, p, Li, Lj, Lk, ms, RTOL, ATOL, NSTEPS)
end;

#= -----------------------------------------------------
Two cases are considered:
Case 1 - Manufactured mesh, `NeoHookean` material and GlobalLoad
Case 2 - GMSH mesh and `HyperElastic` material and Pressure
------------------------------------------------------=#
abstract type AbstractCase end
struct FirstCase <: AbstractCase end
struct SecondCase <: AbstractCase end

"Return the boundary conditions"
function boundary_conditions()
    (; p) = parameters()

    bc_fixed_x_label = "fixed-ux"
    bc_fixed_x = FixedDof(:u, [1], bc_fixed_x_label)
    bc_fixed_y_label = "fixed-uj"
    bc_fixed_y = FixedDof(:u, [2], bc_fixed_y_label)
    bc_fixed_k_label = "fixed-uk"
    bc_fixed_k = FixedDof(:u, [3], bc_fixed_k_label)
    bc_load_label = "compression"
    bc_load = GlobalLoad(:u, t -> [-p * t, 0, 0], bc_load_label)
    bc_labels = [bc_fixed_x_label, bc_fixed_y_label, bc_fixed_k_label, bc_load_label]

    (; bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_load, bc_labels)
end;

"Return the first case structural model"
function structure(::FirstCase=FirstCase())
    (; μ, K, Li, Lj, Lk) = parameters()
    # -------------
    # Mesh
    # -------------
    n₁ = Node(0.0, 0.0, 0.0)
    n₂ = Node(0.0, 0.0, Lk)
    n₃ = Node(0.0, Lj, Lk)
    n₄ = Node(0.0, Lj, 0.0)
    n₅ = Node(Li, 0.0, 0.0)
    n₆ = Node(Li, 0.0, Lk)
    n₇ = Node(Li, Lj, Lk)
    n₈ = Node(Li, Lj, 0.0)
    vec_nodes = [n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈]
    mesh = Mesh(; nodes=vec_nodes)
    f₁ = TriangularFace(n₅, n₈, n₆, "loaded_face_1")
    f₂ = TriangularFace(n₆, n₈, n₇, "loaded_face_2")
    f₃ = TriangularFace(n₄, n₁, n₂, "x=0_face_1")
    f₄ = TriangularFace(n₄, n₂, n₃, "x=0_face_2")
    f₅ = TriangularFace(n₆, n₂, n₁, "y=0_face_1")
    f₆ = TriangularFace(n₆, n₁, n₅, "y=0_face_2")
    f₇ = TriangularFace(n₁, n₄, n₅, "z=0_face_1")
    f₈ = TriangularFace(n₄, n₈, n₅, "z=0_face_2")
    vec_faces = [f₁, f₂, f₃, f₄, f₅, f₆, f₇, f₈]
    append!(faces(mesh), vec_faces)
    ## Entities
    t₁ = Tetrahedron(n₁, n₄, n₂, n₆, "tetra_1")
    t₂ = Tetrahedron(n₆, n₂, n₃, n₄, "tetra_2")
    t₃ = Tetrahedron(n₄, n₃, n₆, n₇, "tetra_3")
    t₄ = Tetrahedron(n₄, n₁, n₅, n₆, "tetra_4")
    t₅ = Tetrahedron(n₄, n₆, n₅, n₈, "tetra_5")
    t₆ = Tetrahedron(n₄, n₇, n₆, n₈, "tetra_6")
    vec_elems = [t₁, t₂, t₃, t₄, t₅, t₆]
    append!(elements(mesh), vec_elems)
    dof_dim = 3
    set_dofs!(mesh, :u, dof_dim)
    # -------------------------------
    # Materials
    # -------------------------------
    neo_hookean = NeoHookean(K, μ, "NeoBuiltIn")
    materials = StructuralMaterial(neo_hookean => [t₁, t₂, t₃, t₄, t₅, t₆])
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    (; bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_load) = boundary_conditions()
    face_bc = [bc_fixed_x => [f₃, f₄],
               bc_fixed_y => [f₅, f₆],
               bc_fixed_k => [f₇, f₈],
               bc_load => [f₁, f₂]]
    bcs = StructuralBoundaryCondition(face_bc)

    Structure(mesh, materials, bcs)
end;

"Return the second case structural model"
function structure(::SecondCase)
    (; μ, K, p, Li, Lj, Lk, ms) = parameters()
    # -------------------------------
    # Materials
    # -------------------------------
    "Neo-Hookean strain energy function given the Green-Lagrange strain
    tensor `𝔼`, second lamé parameter `μ` and bulk modulus `K`."
    function strain_energy_neo(𝔼::AbstractMatrix, K::Real, μ::Real)
        # Right hand Cauchy strain tensor
        ℂ = Symmetric(2 * 𝔼 + eye(3))
        J = sqrt(det(ℂ))
        # First invariant
        I₁ = tr(ℂ)
        # Strain energy function
        Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
    end
    # The order must be the same defined in the strain energy (splatting)
    params = [K, μ]
    mat_label = "neoHyper"
    neo_hookean_hyper = HyperElastic(params, strain_energy_neo, mat_label)
    materials = StructuralMaterial(neo_hookean_hyper)
    # -------------------------------
    # Boundary Conditions
    # -------------------------------
    (; bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_labels) = boundary_conditions()
    bc_load = Pressure(:u, t -> p * t, last(bc_labels))
    bcs = StructuralBoundaryCondition(bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_load)
    # -------------------------------
    # Entities
    # -------------------------------
    faces_label = "triangle"
    elems_label = "tetrahedron"
    vfaces = [TriangularFace(faces_label)]
    velems = [Tetrahedron(elems_label)]
    s_entities = StructuralEntity(velems, vfaces)
    entities_labels = [faces_label, elems_label]
    # -------------------------------
    # Mesh
    # -------------------------------
    filename = "uniaxial_compression"
    labels = [mat_label, entities_labels, bc_labels]
    local mesh_path
    output = @capture_out begin
        mesh_path = create_uniaxial_mesh(Li, Lj, Lk, labels, filename, ms)
    end
    gmsh_println(output)
    msh_file = MshFile(mesh_path)
    mesh = Mesh(msh_file, s_entities)
    dof_dim = 3
    set_dofs!(mesh, :u, dof_dim)
    apply!(materials, mesh)
    apply!(bcs, mesh)

    Structure(mesh, materials, bcs)
end;

"Return the problem solution"
function solve(case::AbstractCase)
    (; NSTEPS) = parameters()
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    s = structure(case)
    sa = NonLinearStaticAnalysis(s; NSTEPS)
    reset!(sa)
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

sa = solve(FirstCase())

#= TODO: REVISIT THIS TESTS

function run_uniaxial_compression()
    "Computes numeric solution α(L_def/L_ref), β(L_def/L_ref) and γ(L_def/L_ref)
    for analytic validation."
    function αβγ_numeric(states_sol::AbstractSolution)
        s = structure(analysis(states_sol))
        # Node at (Li, Lj, Lk)
        n₇ = nodes(s)[7]
        displacements_n₇ = displacements(states_sol_case₁, n₇)
        # Displacements in the x (component 1) axis at node 7
        numerical_uᵢ = displacements_n₇[1]
        numerical_α = 1 .+ numerical_uᵢ / Li
        # Displacements in the y (component 2) axis at node 7
        numerical_uⱼ = displacements_n₇[2]
        numerical_β = 1 .+ numerical_uⱼ / Lj
        # Displacements in the z (component 3) axis at node 7
        numerical_uₖ = displacements_n₇[3]
        numerical_γ = 1 .+ numerical_uₖ / Lk
        numerical_α, numerical_β, numerical_γ, numerical_uᵢ, numerical_uⱼ, numerical_uₖ
    end
    # Numeric solution for testing
    numeric_α_case₁, numeric_β_case₁, numeric_γ_case₁, numeric_uᵢ_case₁, _, _ = αβγ_numeric(states_sol_case₁)
    # Extract ℙ and ℂ from the last state using a random element
    e = rand(elements(s₁))
    # Cosserat or second Piola-Kirchhoff stress tensor
    ℙ_numeric_case₁ = stress(states_sol_case₁, e)
    # ℙᵢᵢ component:
    ℙᵢᵢ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 1, 1)
    # ℙⱼⱼ component:
    ℙⱼⱼ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 2, 2)
    # ℙₖₖ component:
    ℙₖₖ_numeric_case₁ = getindex.(ℙ_numeric_case₁, 3, 3)
    # Get the Right hand Cauchy strain tensor ℂ at a random state
    ℂ_rand_numeric_case₁ = rand(strain(states_sol_case₁, e))
    # Get the Second Piola Kirchhoff stress tensor ℙ at a random state
    ℙ_rand_numeric_case₁ = rand(stress(states_sol_case₁, e))
    # Load factors
    load_factors_case₁ = load_factors(sa)
    # -----------------------------------------------
    # Case 2 - GMSH mesh and `HyperElastic` material
    #------------------------------------------------
    # -------------------------------
    # Materials
    # -------------------------------
    # Define a new HyperElastic material from the strain energy function
    "Neo-Hookean strain energy function given the Green-Lagrange strain
    tensor `𝔼`, second lamé parameter `μ` and bulk modulus `K`."
    function strain_energy_neo(𝔼::AbstractMatrix, K::Real, μ::Real)
        # Right hand Cauchy strain tensor
        ℂ = Symmetric(2 * 𝔼 + eye(3))
        J = sqrt(det(ℂ))
        # First invariant
        I₁ = tr(ℂ)
        # Strain energy function
        Ψ = μ / 2 * (I₁ - 2 * log(J)) + K / 2 * (J - 1)^2
    end
    params = [K, μ] # The order must be the same defined in the strain energy (splatting)
    mat_label = "neoHyper"
    neo_hookean_hyper = HyperElastic(params, strain_energy_neo, mat_label)
    # Material types without assigned elements
    materials = StructuralMaterial(neo_hookean_hyper)
    # -------------------------------
    # Boundary Conditions
    # -------------------------------
    # Redefine the load boundary condition
    bc_load = Pressure(:u, t -> p * t, bc₄_label)
    # BoundaryConditions types without assigned node, feces and elements
    s_boundary_conditions = StructuralBoundaryCondition(bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_load)
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
    filename = "uniaxial_compression"
    labels = [mat_label, entities_labels, bc_labels]
    local mesh_path
    output = @capture_out begin
        mesh_path = create_uniaxial_mesh(Li, Lj, Lk, labels, filename, ms)
    end
    gmsh_println(output)
    msh_file = MshFile(mesh_path)
    mesh = Mesh(msh_file, s_entities)
    # -------------------------------
    # Dofs
    #--------------------------------
    dof_dim = 3
    set_dofs!(mesh, :u, dof_dim)
    # -------------------------------
    # Structure
    # -------------------------------
    apply!(materials, mesh)
    apply!(s_boundary_conditions, mesh)
    s₂ = Structure(mesh, materials, s_boundary_conditions)
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    sa₂ = NonLinearStaticAnalysis(s₂, load_factors(sa))
    reset!(sa₂)
    # -------------------------------
    # Numerical solution
    # -------------------------------
    states_sol_case₂ = solve!(sa₂, nr)
    # Extract ℙ and ℂ from the last state using a random element
    e = rand(elements(s₂))
    # Numeric solution for testing
    numeric_α_case₂, numeric_β_case₂, numeric_γ_case₂, numeric_uᵢ_case₂, _, _ = αβγ_numeric(states_sol_case₂)
    # Cosserat or second Piola-Kirchhoff stress tensor
    ℙ_numeric_case₂ = stress(states_sol_case₂, e)
    # ℙᵢᵢ component:
    ℙᵢᵢ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 1, 1)
    # ℙⱼⱼ component:
    ℙⱼⱼ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 2, 2)
    # ℙₖₖ component:
    ℙₖₖ_numeric_case₂ = getindex.(ℙ_numeric_case₂, 3, 3)
    # Get the Right hand Cauchy strain tensor ℂ at a random state
    ℂ_rand_numeric_case₂ = rand(strain(states_sol_case₂, e))
    # Get the Second Piola Kirchhoff stress tensor ℙ at a random state
    ℙ_rand_numeric_case₂ = rand(stress(states_sol_case₂, e))
    # Load factors
    load_factors_case₂ = load_factors(sa₂)
    #-----------------------------
    # Analytic solution
    #-----------------------------
    "Computes displacements numeric solution uᵢ, uⱼ and uₖ for analytic validation."
    function u_ijk_numeric(numerical_α::Vector{<:Real}, numerical_β::Vector{<:Real},
                           numerical_γ::Vector{<:Real},
                           x::Real, y::Real, z::Real)
        x * (numerical_α .- 1), y * (numerical_β .- 1), z * (numerical_γ .- 1)
    end
    # Test with Second Piola-Kirchoff stress tensor `ℙ`.
    "Computes ℙ(1,1) given α, β and γ."
    function analytic_ℙᵢᵢ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        μ * α - μ * (α .^ (-1)) + K * (β .^ 2) .* (α .* (β .^ 2) .- 1)
    end
    "Computes ℙ(2,2) given α, β and γ."
    function analytic_ℙⱼⱼ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        μ * β - μ * (β .^ (-1)) + K * β .* ((α .^ 2) .* (β .^ 2) - α)
    end
    "Computes ℙ(2,2) given α, β and γ."
    function analytic_ℙₖₖ(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        analytic_ℙⱼⱼ(α, β, μ, K)
    end
    # Compute the analytic Second Piola-Kirchoff stress tensor `ℙ` for the numeric vectors α and β
    # Case 1
    ℙᵢᵢ_analytic_case₁ = analytic_ℙᵢᵢ(numeric_α_case₁, numeric_β_case₁)
    ℙⱼⱼ_analytic_case₁ = analytic_ℙⱼⱼ(numeric_α_case₁, numeric_β_case₁)
    ℙₖₖ_analytic_case₁ = analytic_ℙₖₖ(numeric_α_case₁, numeric_β_case₁)
    # Case 2
    ℙᵢᵢ_analytic_case₂ = analytic_ℙᵢᵢ(numeric_α_case₂, numeric_β_case₂)
    ℙⱼⱼ_analytic_case₂ = analytic_ℙⱼⱼ(numeric_α_case₂, numeric_β_case₂)
    ℙₖₖ_analytic_case₂ = analytic_ℙₖₖ(numeric_α_case₂, numeric_β_case₂)
    # -------------------------------
    # Interpolator tests for Case 2
    #--------------------------------
    rand_point = [[rand() * Li, rand() * Lj, rand() * Lk]]
    eval_handler_rand = PointEvalHandler(mesh(s₂), rand_point)
    # Compute analytic solution at a random point
    uᵢ_case₂, uⱼ_case₂, uₖ_case₂ = u_ijk_numeric(numeric_α_case₂, numeric_β_case₂, numeric_γ_case₂,
                                                 rand_point[]...)
    rand_point_uᵢ = displacements(states_sol_case₂, eval_handler_rand, 1)
    rand_point_uⱼ = displacements(states_sol_case₂, eval_handler_rand, 2)
    rand_point_uₖ = displacements(states_sol_case₂, eval_handler_rand, 3)
    stress_point = stress(states_sol_case₂, eval_handler_rand)[]
    #-----------------------------
    # Test boolean for CI
    #-----------------------------
    @testset "Case 1 Uniaxial Compression Example" begin
        @test ℙᵢᵢ_analytic_case₁ ≈ ℙᵢᵢ_numeric_case₁ rtol = RTOL
        @test ℙᵢᵢ_analytic_case₁ ≈ ℙᵢᵢ_numeric_case₁ rtol = RTOL
        @test ℙⱼⱼ_analytic_case₁ ≈ ℙⱼⱼ_numeric_case₁ atol = ATOL
        @test ℙₖₖ_analytic_case₁ ≈ ℙₖₖ_numeric_case₁ atol = ATOL
        @test norm(ℙⱼⱼ_analytic_case₁) ≈ 0 atol = ATOL
        @test norm(ℙₖₖ_analytic_case₁) ≈ 0 atol = ATOL
        @test -p * load_factors_case₁ ≈ ℙᵢᵢ_analytic_case₁ rtol = RTOL
    end

    @testset "Case 2 Uniaxial Compression Example" begin
        @test ℙᵢᵢ_analytic_case₂ ≈ ℙᵢᵢ_numeric_case₂ rtol = RTOL
        @test ℙᵢᵢ_analytic_case₂ ≈ ℙᵢᵢ_numeric_case₂ rtol = RTOL
        @test ℙⱼⱼ_analytic_case₂ ≈ ℙⱼⱼ_numeric_case₂ atol = ATOL
        @test ℙₖₖ_analytic_case₂ ≈ ℙₖₖ_numeric_case₂ atol = ATOL
        @test norm(ℙⱼⱼ_analytic_case₂) ≈ 0 atol = ATOL
        @test norm(ℙₖₖ_analytic_case₂) ≈ 0 atol = ATOL
        @test -p * load_factors_case₂ ≈ ℙᵢᵢ_analytic_case₂ rtol = RTOL
        # Interpolation
        @test uᵢ_case₂ ≈ rand_point_uᵢ rtol = RTOL
        @test uⱼ_case₂ ≈ rand_point_uⱼ rtol = RTOL
        @test uₖ_case₂ ≈ rand_point_uₖ rtol = RTOL
        @test getindex.(stress_point, 1) ≈ ℙᵢᵢ_analytic_case₂ rtol = RTOL
    end
end

run_uniaxial_compression()

=#
