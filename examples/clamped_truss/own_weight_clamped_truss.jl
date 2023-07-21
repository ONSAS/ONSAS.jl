# ---------------------------------
# Clamped truss with self-weight
# --------------------------------
#=
This model consist of a clamped truss with its own weight. First the
weight is applied and consecutive the load at the time starting from the own
weight deformed configuration.

# Analysis 1: Own weight

## Forces
-> g
--------> F

# Analysis 1: Own weight
#-----------------------------
# ->  ->  ->  ->  ->  ->  -> -
#-----------------------------

# Analysis 2: Own weight + load
#-----------------------------
# ->  ->  ->  ->  ->  ->  -> - -----> F
#-----------------------------
=#
using ONSAS, Test
using Dictionaries: dictionary

"Return problem parameters"
function parameters()
    N = 10     # Number of elements.
    E = 30e6    # Young's modulus.
    ν = 0.3     # Poisson's ratio.
    ρ = 750.0   # Density.
    L = 200     # Element length.
    A = 1       # Cross section area.
    F = 10e6    # Force at the tip
    g = 9.81    # Gravity
    ϵ_model = RotatedEngineeringStrain
    g, N, E, ν, ρ, L, A, F, ϵ_model
end

#-----------------------------
# Analytic solution
#-----------------------------
"Analytic displacement uᵢ at the tip"
function analytic_u(x::Real; F::Real, E::Real,
                    A::Real, L::Real, g::Real, ρ::Real)
    b = ρ * g
    C = F / (E * A) + b * L / (E * A)
    -b / (2 * E * A) * x^2 + C * x
end

"Return mesh"
function mesh(N::Int, L::Real, A::Real, ϵ_model)
    nodes = [Node(l) for l in LinRange(0, L, N + 1)]
    elements = [Truss(nodes[i], nodes[i + 1], Square(sqrt(A)), ϵ_model)
                for i in 1:N]
    mesh = Mesh(; nodes, elements)
    dof_dim = 1
    set_dofs!(mesh, :u, dof_dim)
    mesh
end;

"Return structure"
function structure(N::Int;
                   E::Real, ν::Real, ρ::Real, L::Real, A::Real, g::Real, ϵ_model)
    # -------------
    # Mesh
    # -------------
    mesh = Main.mesh(N, L, A, ϵ_model)
    # -------------------------------
    # Materials
    # -------------------------------
    material = SVK(; E, ν, ρ, label="material")
    materials = StructuralMaterial(material => elements(mesh))
    # -------------------------------
    # Boundary conditions
    # -------------------------------
    fixed_bc = FixedDof(:u, [1], "fixed")
    gravity_bc_ramp = GlobalLoad(:u, t -> t * [density(material) * g], "gravity")
    node_bcs = dictionary([fixed_bc => [first(nodes(mesh))]])
    element_bcs = dictionary([gravity_bc_ramp => elements(mesh)])
    boundary_conditions = StructuralBoundaryCondition(; node_bcs, element_bcs)
    # -------------------------------
    # Structure
    # -------------------------------
    Structure(mesh, materials, boundary_conditions)
end;
#-----------------------------
# Problem parameters
#-----------------------------
g, N, E, ν, ρ, L, A, F, ϵ_model = parameters()
#-----------------------------
# Analysis 1: Own weight
#-----------------------------
# Structure
s = structure(N; E, ν, ρ, L, A, g, ϵ_model)
# Analysis
gravity_analysis = LinearStaticAnalysis(s; NSTEPS=10)
# Solution
gravity_solution = solve!(gravity_analysis)
# Tests
# External force
numerical_Fext = external_forces(last(states(gravity_solution)))
analytical_Fext = [ρ * g * A * L / N for n in nodes(s)]
analytical_Fext[1] = analytical_Fext[end] = ρ * g * A / 2 * L / N
@test numerical_Fext ≈ analytical_Fext atol = 1e-6
# Displacement
numerical_u_last_node = last(displacements(gravity_solution, last(nodes(s)), 1))
analytic_u_last_node = analytic_u(L; F=0.0, E, A, L, g, ρ)
@test numerical_u_last_node ≈ analytic_u_last_node atol = 1e-6
#-----------------------------------
# Analysis 2: Own weight + tip load
#-----------------------------------
# tip_load_bc = GlobalLoad(:u, t -> t * [F], "tip_load")
constant_gravity = GlobalLoad(:u, t -> [ρ * g], "gravity")
replace!(boundary_conditions(s), constant_gravity)
