# -----------------------------------------
# Hyperelastic Uniaxial Compression Example
# -----------------------------------------
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
    ms = 0.5                    # Mesh size parameter
    RTOL = 1e-4                 # Relative tolerance for tests
    ATOL = 1e-8                # Absolute tolerance for tests
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
    bc_fixed_x = FixedField(:u, [1], bc_fixed_x_label)
    bc_fixed_y_label = "fixed-uj"
    bc_fixed_y = FixedField(:u, [2], bc_fixed_y_label)
    bc_fixed_k_label = "fixed-uk"
    bc_fixed_k = FixedField(:u, [3], bc_fixed_k_label)
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
    n1 = Node(0.0, 0.0, 0.0)
    n2 = Node(0.0, 0.0, Lk)
    n3 = Node(0.0, Lj, Lk)
    n4 = Node(0.0, Lj, 0.0)
    n5 = Node(Li, 0.0, 0.0)
    n6 = Node(Li, 0.0, Lk)
    n7 = Node(Li, Lj, Lk)
    n8 = Node(Li, Lj, 0.0)
    vec_nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
    mesh = Mesh(; nodes=vec_nodes)
    f1 = TriangularFace(n5, n8, n6, "loaded_face_1")
    f2 = TriangularFace(n6, n8, n7, "loaded_face_2")
    f3 = TriangularFace(n4, n1, n2, "x=0_face_1")
    f4 = TriangularFace(n4, n2, n3, "x=0_face_2")
    f5 = TriangularFace(n6, n2, n1, "y=0_face_1")
    f6 = TriangularFace(n6, n1, n5, "y=0_face_2")
    f7 = TriangularFace(n1, n4, n5, "z=0_face_1")
    f8 = TriangularFace(n4, n8, n5, "z=0_face_2")
    vec_faces = [f1, f2, f3, f4, f5, f6, f7, f8]
    append!(faces(mesh), vec_faces)
    ## Entities
    t1 = Tetrahedron(n1, n4, n2, n6, "tetra_1")
    t2 = Tetrahedron(n6, n2, n3, n4, "tetra_2")
    t3 = Tetrahedron(n4, n3, n6, n7, "tetra_3")
    t4 = Tetrahedron(n4, n1, n5, n6, "tetra_4")
    t5 = Tetrahedron(n4, n6, n5, n8, "tetra_5")
    t6 = Tetrahedron(n4, n7, n6, n8, "tetra_6")
    vec_elems = [t1, t2, t3, t4, t5, t6]
    append!(elements(mesh), vec_elems)
    dof_dim = 3
    set_dofs!(mesh, :u, dof_dim)
    # -------------------------------
    # Materials
    # -------------------------------
    neo_hookean = NeoHookean(K, μ, "NeoBuiltIn")
    materials = StructuralMaterial(neo_hookean => [t1, t2, t3, t4, t5, t6])
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    (; bc_fixed_x, bc_fixed_y, bc_fixed_k, bc_load) = boundary_conditions()
    face_bc = [bc_fixed_x => [f3, f4],
               bc_fixed_y => [f5, f6],
               bc_fixed_k => [f7, f8],
               bc_load => [f1, f2]]
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
    ONSAS.solve(sa, nr)
end;

"Computes numeric solution
α(L_def/L_ref), β(L_def/L_ref) and γ(L_def/L_ref)
for analytic validation."
function αβγ_numeric(sol::AbstractSolution)
    (; Li, Lj, Lk) = parameters()
    s = ONSAS.structure(analysis(sol))
    # Node at (Li, Lj, Lk)
    n7 = nodes(s)[7]
    displacements_n7 = displacements(sol, n7)
    # Displacements in the x (component 1) axis at node 7
    ux = displacements_n7[1]
    α = 1 .+ ux / Li
    # Displacements in the y (component 2) axis at node 7
    uy = displacements_n7[2]
    β = 1 .+ uy / Lj
    # Displacements in the z (component 3) axis at node 7
    uz = displacements_n7[3]
    γ = 1 .+ uz / Lk
    α, β, γ, ux, uy, uz
end;

"Return relevant numeric results for testing"
function numerical_solution(sol::AbstractSolution)
    s = ONSAS.structure(analysis(sol))
    # Numeric solution for testing
    α, β, γ, uᵢ, _, _ = αβγ_numeric(sol)
    # Extract P and ℂ from the last state using a random element
    e = rand(elements(s))
    # Cosserat or second Piola-Kirchhoff stress tensor
    P = stress(sol, e)
    # P11 component:
    P11 = getindex.(P, 1, 1)
    # P22 component:
    P22 = getindex.(P, 2, 2)
    # P33 component:
    P33 = getindex.(P, 3, 3)
    # Get the Right hand Cauchy strain tensor C at a random state
    C_rand = rand(strain(sol, e))

    P11, P22, P33, α, β, γ, uᵢ, P
end;

"Computes displacements numeric solution uᵢ, uⱼ and uₖ for analytic validation."
function u_ijk_numeric(numerical_α::Vector{<:Real},
                       numerical_β::Vector{<:Real},
                       numerical_γ::Vector{<:Real},
                       x::Real, y::Real, z::Real)
    x * (numerical_α .- 1), y * (numerical_β .- 1), z * (numerical_γ .- 1)
end;

"Return analytic Piola-Kirchoff stress tensor `P`"
function analytic_P(sol::AbstractSolution)
    _, _, _, numeric_α, numeric_β, _, _ = numerical_solution(sol::AbstractSolution)
    (; μ, K) = parameters()

    # Test with Second Piola-Kirchoff stress tensor `P`.
    "Computes P(1,1) given α, β and γ."
    function analytic_P11(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        μ * α - μ * (α .^ (-1)) + K * (β .^ 2) .* (α .* (β .^ 2) .- 1)
    end
    "Computes P(2,2) given α, β and γ."
    function analytic_P22(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        μ * β - μ * (β .^ (-1)) + K * β .* ((α .^ 2) .* (β .^ 2) - α)
    end
    "Computes P(2,2) given α, β and γ."
    function analytic_P33(α::Vector{<:Real}, β::Vector{<:Real}, μ::Real=μ, K::Real=K)
        analytic_P22(α, β, μ, K)
    end
    P11_analytic = analytic_P11(numeric_α, numeric_β)
    P22_analytic = analytic_P22(numeric_α, numeric_β)
    P33_analytic = analytic_P33(numeric_α, numeric_β)

    P11_analytic, P22_analytic, P33_analytic
end;

"Test case numerical example solution"
function test(sol::AbstractSolution)
    (; RTOL, ATOL, p, Li, Lj, Lk) = parameters()
    a = analysis(sol)
    P11_numeric, P22_numeric, P33_numeric, _, _ = numerical_solution(sol)
    P11_analytic, P22_analytic, P33_analytic = analytic_P(sol)
    @test P11_analytic ≈ P11_numeric rtol = RTOL
    @test P22_analytic ≈ P22_numeric atol = ATOL
    @test P33_analytic ≈ P33_numeric atol = ATOL
    @test -p * load_factors(a) ≈ P11_analytic rtol = RTOL
    @test norm(P22_numeric) ≈ 0 atol = ATOL
    @test norm(P33_numeric) ≈ 0 atol = ATOL

    rand_point = [rand() * [Li, Lj, Lk]]
    ph = PointEvalHandler(mesh(ONSAS.structure(a)), rand_point)

    # Compute analytic solution at a random point
    _, _, _, numeric_α, numeric_β, numeric_γ, _ = numerical_solution(sol)

    u1_case, u2_case, u3_case = u_ijk_numeric(numeric_α, numeric_β, numeric_γ,
                                              rand_point[]...)
    rand_point_u1 = displacements(sol, ph, 1)
    rand_point_u2 = displacements(sol, ph, 2)
    rand_point_u3 = displacements(sol, ph, 3)
    stress_point = stress(sol, ph)[]
    @test u1_case ≈ rand_point_u1 rtol = RTOL
    @test u2_case ≈ rand_point_u2 rtol = RTOL
    @test u3_case ≈ rand_point_u3 rtol = RTOL
    @test getindex.(stress_point, 1) ≈ P11_analytic rtol = RTOL
end;

"Run the example"
function run()
    for case in (FirstCase(), SecondCase())
        sol = solve(case)
        test(sol)
    end
end;

run()
