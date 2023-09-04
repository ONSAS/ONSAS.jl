# ---------------------
# Clamped truss example
# ---------------------
#=
This model is a static generalization taken from [3].
[3] https://github.com/JuliaReach/SetPropagation-FEM-Examples/blob/main/examples/Clamped/Clamped_Model.jl
=#
using Test, LinearAlgebra
using ONSAS

"Return the problem parameters"
function parameters()
    N = 100               # Number of elements.
    E = 30e6              # Young's modulus.
    ν = 0.3               # Poisson's ratio.
    ρ = 7.3e-4            # Density.
    L = 200               # Element length.
    A = 1                 # Cross section area.
    F = 10e6              # Force at the tip
    ϵ_model = GreenStrain # Strain model
    NSTEPS = 10           # Number of load factors steps

    (; NSTEPS, ϵ_model, N, E, ν, ρ, L, A, F)
end;

"Return the problem structural model"
function structure()
    (; N, E, ν, ρ, L, A, F, ϵ_model) = parameters()
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
    steel = SVK(; E=E, ν=ν, ρ=ρ, label="steel")
    materials = StructuralMaterial(steel => elements)
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    bc₁ = FixedDof(:u, [1], "fixed_uₓ")
    bc₂ = GlobalLoad(:u, t -> [F * t], "load in j")
    # Apply bcs to the nodes
    boundary_conditions = StructuralBoundaryCondition(bc₁ => [first(nodes)], bc₂ => [last(nodes)])

    Structure(mesh, materials, boundary_conditions)
end;

"Return the problem solution"
function solve()
    s = structure()
    # -------------------------------
    # Structural Analysis
    # -------------------------------
    (; NSTEPS) = parameters()
    sa = NonLinearStaticAnalysis(s; NSTEPS=NSTEPS)
    # -------------------------------
    # Numerical solution
    # -------------------------------
    solve!(sa, NewtonRaphson())
end;

"Test problem solution"
function test(sol::AbstractSolution)
    (; F, ϵ_model, E, A, L) = parameters()
    # Force and displacement at the tip
    sa = analysis(sol)
    vec_nodes = ONSAS.nodes(mesh(ONSAS.structure(sa)))
    numeric_uᵢ = displacements(sol, last(vec_nodes))[1]
    numeric_P_tip = F * load_factors(sa)
    #-----------------------------
    # Analytic solution
    #-----------------------------
    # Compute the analytic values for the strain, stress and force at the tip
    "Analytic force given `uᵢ` towards x axis at the tip node"
    function analytic_P(::Type{GreenStrain}, uᵢ::Real, E::Real=E, l₀::Real=L, A₀::Real=A)
        ϵ_green = 0.5 * ((l₀ + uᵢ)^2 - l₀^2) / (l₀^2)
        # Cosserat stress
        𝐒₁₁ = E * ϵ_green
        # Piola stress
        𝐏₁₁ = (l₀ + uᵢ) / l₀ * 𝐒₁₁
        𝐏₁₁ * A₀
    end
    #
    analytic_P_tip = analytic_P.(Ref(ϵ_model), numeric_uᵢ)
    @testset "Piola-Kirchoff tensor at the right-most node" begin
        @test analytic_P_tip ≈ numeric_P_tip rtol = 1e-3
    end
end

"Run the example"
function run()
    sol = solve()
    test(sol)
end

run()
