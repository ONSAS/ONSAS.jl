# --------------------- 
# Clamped truss example
# ---------------------
#=
This model is taken from [1]. The theoretical derivations of the analytic solution can be found in [2].
The implementation of stiffness and mass matrices in Julia can be found in [3].

[1] Malakiyeh, Mohammad Mahdi, Saeed Shojaee, and Klaus-Jürgen Bathe. "The Bathe time integration method revisited for prescribing desired numerical dissipation." Computers & Structures 212 (2019): 289-298.

[2] Mechanical Vibrations, Gerardin et al, page 250-251.

[3] https://github.com/JuliaReach/SetPropagation-FEM-Examples/blob/main/examples/Clamped/Clamped_Model.jl
=#
using Test, LinearAlgebra, Dictionaries
using ONSAS

"Runs the clamped truss example."
function run_clamped_truss_example()
    # Parameters
    N = 100     # Number of elements.
    E = 30e6    # Young's modulus.
    ν = 0.3     # Poisson's ratio. 
    ρ = 7.3e-4  # Density.
    L = 200     # Element length.
    A = 1       # Cross section area.
    F = 10e6    # Force at the tip
    ϵ_model = rand((GreenStrain, RotatedEngineeringStrain))
    # -------------
    # Mesh
    # -------------
    nodes = [Node(l) for l in LinRange(0, L, N + 1)]
    elements = [Truss(nodes[i], nodes[i + 1], Square(sqrt(A)), ϵ_model)
                for i in 1:N]
    mesh = Mesh(; nodes, elements)
    # -------------------------------
    # Dofs
    #--------------------------------
    dof_dim = 1
    apply!(mesh, :u, dof_dim)
    # -------------------------------
    # Materials
    # -------------------------------
    steel = Svk(; E=E, ν=ν, ρ=ρ, label="steel")
    mat_dict = dictionary([steel => elements])
    materials = StructuralMaterial(mat_dict)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    # Fixed dofs
    bc₁ = FixedDof(; components=[1], name="fixed_uₓ")
    # Load 
    bc₂ = GlobalLoad(; values=t -> [F * t], name="load in j")
    # Apply bcs to the nodes
    node_bc = dictionary([bc₁ => [first(nodes)], bc₂ => [last(nodes)]])
    boundary_conditions = StructuralBoundaryCondition(; node_bcs=node_bc)
    # -------------------------------
    # Structure
    # -------------------------------
    s = Structure(mesh, materials, boundary_conditions)
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    NSTEPS = 10
    sa = NonLinearStaticAnalysis(s; NSTEPS=NSTEPS)
    # -------------------------------
    # Algorithm
    # -------------------------------
    # nr default tolerances
    alg = NewtonRaphson()
    # -------------------------------
    # Numerical solution
    # -------------------------------
    solution = solve!(sa, alg)
    # Extract displacement at the tip
    numeric_uᵢ = displacements(solution, last(nodes))[1]
    numeric_σ_tip_tensor = stress(solution, last(elements))
    numeric_σ_tip = getindex.(numeric_σ_tip_tensor, 1)
    numeric_ϵ_tip_tensor = strain(solution, last(elements))
    numeric_ϵ_tip = getindex.(numeric_ϵ_tip_tensor, 1)
    numeric_F_tip = F * load_factors(sa)
    #-----------------------------
    # Analytic solution  
    #-----------------------------
    # Compute the analytic values for the strain, stress and force at the tip
    "Analytic rotated engineering strain solution for the displacement 
    `uᵢ` towards x axis at the tip node."
    function analytic_ϵ(::Type{RotatedEngineeringStrain}, uᵢ::Real, l₀::Real=L)
        ((l₀ + uᵢ)^2 - l₀^2) / (l₀ * (l₀ + (l₀ + uᵢ)))
    end
    function analytic_ϵ(::Type{GreenStrain}, uᵢ::Real, l₀::Real=L)
        ((l₀ + uᵢ)^2 - l₀^2) / (2 * l₀^2)
    end
    "Analytic stress value for a given strain `ϵ`."
    analytic_σ(analytic_ϵ::Vector{<:Real}, E::Real=E) = analytic_ϵ * E
    "Analytic force value for a given strain `ϵ`."
    analytic_F(analytic_σ::Vector{<:Real}, A₀::Real=A) = analytic_σ * A₀
    #
    analytic_ϵ_tip = analytic_ϵ.(ϵ_model, numeric_uᵢ)
    analytic_σ_tip = analytic_σ(analytic_ϵ_tip, E)
    analytic_F_tip = analytic_F(analytic_σ_tip, A)
    #-----------------------------
    # Test boolean for CI  
    #-----------------------------
    @test analytic_F_tip ≈ numeric_F_tip rtol = 1e-3
    @test numeric_σ_tip ≈ analytic_σ_tip rtol = 1e-3
    @test numeric_ϵ_tip ≈ analytic_ϵ_tip rtol = 1e-3
end

run_clamped_truss_example()
