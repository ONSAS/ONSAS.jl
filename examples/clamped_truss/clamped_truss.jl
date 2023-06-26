# ---------------------
# Clamped truss example
# ---------------------
#=
This model is a static generalization taken from [3].
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
    ϵ_model = GreenStrain
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
    set_dofs!(mesh, :u, dof_dim)
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
    bc₁ = FixedDof(:u, [1], "fixed_uₓ")
    # Load
    bc₂ = GlobalLoad(:u, t -> [F * t], "load in j")
    # Apply bcs to the nodes
    boundary_conditions = StructuralBoundaryCondition(bc₁ => [first(nodes)], bc₂ => [last(nodes)])
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
    # Solver
    # -------------------------------
    nr = NewtonRaphson()
    # -------------------------------
    # Numerical solution
    # -------------------------------
    solution = solve!(sa, nr)
    # Force and displacement at the tip
    numeric_uᵢ = displacements(solution, last(nodes))[1]
    numeric_F_tip = F * load_factors(sa)
    #-----------------------------
    # Analytic solution
    #-----------------------------
    # Compute the analytic values for the strain, stress and force at the tip
    "Analytic force given `uᵢ` towards x axis at the tip node."
    function analytic_F(::Type{GreenStrain}, uᵢ::Real, E::Real=E, l₀::Real=L, A₀::Real=A)
        ϵ_green = 0.5 * ((l₀ + uᵢ)^2 - l₀^2) / (l₀^2)
        # Cosserat stress
        𝐒₁₁ = E * ϵ_green
        # Piola stress
        𝐏₁₁ = (l₀ + uᵢ) / l₀ * 𝐒₁₁
        𝐏₁₁ * A₀
    end
    #
    analytic_F_tip = analytic_F.(Ref(ϵ_model), numeric_uᵢ)
    #-----------------------------
    # Test boolean for CI
    #-----------------------------
    @test analytic_F_tip ≈ numeric_F_tip rtol = 1e-3
end

run_clamped_truss_example()
