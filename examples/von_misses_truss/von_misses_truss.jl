# ------------------------------------------------------------- 
# Von Misses Truss Example from (Zerpa, Bazzano 2017 ) - 2.5.4
# -------------------------------------------------------------
using Test, LinearAlgebra
using ONSAS

"Runs the Von Misses Truss example."
function run_von_misses_truss_example()
    ## scalar parameters
    E = 210e9                 # Young modulus in Pa
    ν = 0.0                   # Poisson's modulus
    A₀ = 2.5e-3                # Cross-section area in m²
    ANG = 65                  # truss angle in degrees
    L = 2                     # Length in m 
    V = L * cos(deg2rad(ANG)) # vertical distance in m 
    H = L * sin(deg2rad(ANG)) # horizontal distance in m
    Fₖ = -3e8                 # Vertical load in N
    RTOL = 1e-4               # Relative tolerance for tests
    # -------------
    # Mesh
    # -------------
    ## Nodes
    n₁ = Node(0.0, 0.0, 0.0)
    n₂ = Node(V, 0.0, H)
    n₃ = Node(2V, 0.0, 0.0)
    vec_nodes = [n₁, n₂, n₃]
    ## Cross sections
    d = sqrt(4 * A₀ / pi)
    s₁ = Circle(d)
    a = sqrt(A₀)
    s₂ = Square(a)
    ## Elements 
    truss₁ = Truss(n₁, n₂, s₁, "left_truss") # [n₁, n₂]
    truss₂ = Truss(n₂, n₃, s₂, "right_truss") # [n₂, n₃]
    vec_elems = [truss₁, truss₂]
    ## Mesh
    s_mesh = Mesh(vec_nodes, vec_elems)
    # -------------------------------
    # Dofs
    #--------------------------------
    dof_dim = 3
    apply!(s_mesh, :u, dof_dim)
    # -------------------------------
    # Materials
    # -------------------------------
    steel = SVK(E=E, ν=ν, label="steel")
    mat_dict = dictionary([steel => [truss₁, truss₂]])
    s_materials = StructuralMaterials(mat_dict)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    # Fixed dofs
    bc₁ = FixedDofBoundaryCondition([:u], [1, 2, 3], "fixed_uₓ_uⱼ_uₖ")
    bc₂ = FixedDofBoundaryCondition([:u], [2], "fixed_uⱼ")
    # Load 
    bc₃ = GlobalLoadBoundaryCondition([:u], t -> [0, 0, Fₖ * t], "load in j")
    node_bc = dictionary([bc₁ => [n₁, n₃], bc₂ => [n₂], bc₃ => [n₂]])
    s_boundary_conditions = StructuralBoundaryConditions(node_bcs=node_bc)
    # -------------------------------
    # Structure
    # -------------------------------
    s = Structure(s_mesh, s_materials, s_boundary_conditions)
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    # Final load factor
    λ₁ = 1
    NSTEPS = 10
    sa = NonLinearStaticAnalysis(s, λ₁, NSTEPS=NSTEPS)
    # -------------------------------
    # Algorithm
    # -------------------------------
    tol_f = 1e-7
    tol_u = 1e-7
    max_iter = 100
    tols = ConvergenceSettings(tol_u, tol_f, max_iter)
    alg = NewtonRaphson(tols)
    # -------------------------------
    # Numerical solution
    # -------------------------------
    states_sol = solve!(sa, alg)
    n₂_displacements = displacements(states_sol, n₂)
    numerical_uᵢ = n₂_displacements[1]
    numerical_uⱼ = n₂_displacements[2]
    numerical_uₖ = n₂_displacements[3]
    @test norm(numerical_uᵢ) ≤ RTOL
    @test norm(numerical_uⱼ) ≤ RTOL
    numerical_λᵥ = -load_factors(sa) * Fₖ
    # Test stress and strains 
    σ_truss₂ = stress(states_sol, truss₂)
    ϵ_truss₂ = strain(states_sol, truss₂)
    @test σ_truss₂ == E * ϵ_truss₂
    #-----------------------------
    # Analytic solution  
    #-----------------------------
    "Analytic load factor solution for the displacement `uₖ` towards z axis at node `n₂`."
    function load_factors_analytic(uₖ::Real, E::Real=E, A::Real=A₀, H::Real=H, V::Real=V, l₀=L)
        λ = -2 * E * A *
            ((H + uₖ)^2 + V^2 - l₀^2) /
            (l₀ * (l₀ + sqrt((H + uₖ)^2 + V^2))) *
            (H + uₖ) / sqrt((H + uₖ)^2 + V^2)
    end
    analytics_λᵥ = load_factors_analytic.(numerical_uₖ)
    #-----------------------------
    # Test boolean for CI  
    #-----------------------------
    @test analytics_λᵥ ≈ numerical_λᵥ rtol = RTOL
end

run_von_misses_truss_example()
