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
using Test: @test
using LinearAlgebra: norm
using ONSAS.StaticAnalyses
# Parameters
N = 1000    # Number of elements.
E = 30e6    # Young's modulus.
ν = 0.3     # Poisson's ratio. 
ρ = 7.3e-4  # Density.
L = 200     # Element length.
A = 1       # Cross section area.
F = 10e3    # Force at the tip
# -------------
# Mesh
# -------------
v_nodes = [Node(l) for l in LinRange(0, L, N + 1)]
v_elements = [Truss(v_nodes[i], v_nodes[i+1], Square(sqrt(A))) for i in 1:N]
s_mesh = Mesh(v_nodes, v_elements)
# -------------------------------
# Dofs
#--------------------------------
dof_dim = 1
add!(s_mesh, :u, dof_dim)
# -------------------------------
# Materials
# -------------------------------
steel = SVK(E, ν, ρ, "steel")
mat_dict = dictionary([steel => v_elements])
s_materials = StructuralMaterials(mat_dict)
# -------------------------------
# Boundary conditions
# -------------------------------
# Fixed dofs
bc₁ = FixedDofBoundaryCondition([:u], [1], "fixed_uₓ")
# Load 
bc₂ = GlobalLoadBoundaryCondition([:u], t -> [F * t], "load in j")
# Apply bcs to the nodes
node_bc = dictionary([bc₁ => [first(v_nodes)], bc₂ => [last(v_nodes)]])
s_boundary_conditions = StructuralBoundaryConditions(node_bcs=node_bc)
# -------------------------------
# Structure
# -------------------------------
s = Structure(s_mesh, s_materials, s_boundary_conditions)
# -------------------------------
# Structural Analysis
# -------------------------------
NSTEPS = 10
sa = StaticAnalysis(s, NSTEPS=NSTEPS)
# -------------------------------
# Algorithm
# -------------------------------
# nr default tolerances
alg = NewtonRaphson()
# -------------------------------
# Numerical solution
# -------------------------------
states_sol = solve(sa, alg)
# Re compute the analysis solution just for fun 
reset!(sa)
states_sol = solve(sa, alg)
# Extract displacement at the tip
numeric_uᵢ = displacements(states_sol, last(v_nodes))[1]
numeric_σ_tip_tensor = stress(states_sol, last(v_elements))
numeric_σ_tip = getindex.(numeric_σ_tip_tensor, 1)
numeric_ϵ_tip_tensor = strain(states_sol, last(v_elements))
numeric_ϵ_tip = getindex.(numeric_ϵ_tip_tensor, 1)
numeric_F_tip = F * load_factors(sa)
#-----------------------------
# Analytic solution  
#-----------------------------
# Compute the analytic values for the strain, stress and force at the tip
"Analytic rotated engineering strain solution for the displacement `uᵢ` towards x axis at the tip node."
analytic_ϵ(uᵢ::Real, l₀::Real=L) = ((l₀ + uᵢ)^2 - l₀^2) / (l₀ * (l₀ + (l₀ + uᵢ)))
"Analytic stress value for a given strain `ϵ`."
analytic_σ(analytic_ϵ::Vector{<:Real}, E::Real=E) = analytic_ϵ * E
"Analytic force value for a given strain `ϵ`."
analytic_F(analytic_σ::Vector{<:Real}, A₀::Real=A) = analytic_σ * A₀
#
analytic_ϵ_tip = analytic_ϵ.(numeric_uᵢ)
analytic_σ_tip = analytic_σ(analytic_ϵ_tip, E)
analytic_F_tip = analytic_F(analytic_σ_tip, A)
#-----------------------------
# Test boolean for CI  
#-----------------------------
@test analytic_F_tip ≈ numeric_F_tip rtol = 1e-3
@test numeric_σ_tip ≈ analytic_σ_tip rtol = 1e-3
@test numeric_ϵ_tip ≈ analytic_ϵ_tip rtol = 1e-3



