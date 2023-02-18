"""
Module defining structural analyses that can be solved. 
Each structural analysis consists of a data type with an structure to be analyzed, analysis parameters and 
an specific state which describes the current state of structure in the analysis.
"""
module StructuralAnalyses

using Reexport: @reexport

@reexport using ..StructuralSolvers
@reexport using ..StructuralModel
@reexport using ..Materials
import ..Elements: internal_forces, inertial_forces, strain, stress
import ..StructuralModel: free_dofs
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
using ..Utils: row_vector

export AbstractStructuralState, assembler, displacements, Δ_displacements,
    residual_forces, iteration_residuals, tangent_matrix, external_forces, structure
export AbstractStructuralAnalysis, initial_time, current_time, final_time,
    _next!, is_done, current_state, current_iteration
export _apply!

""" Abstract supertype to define a new structural state.
**Common methods:**
* [`assembler`](@ref)
* [`displacements`](@ref)
* [`Δ_displacements`](@ref)
* [`external_forces`](@ref)
* [`iteration_residuals`](@ref)
* [`internal_forces`](@ref)
* [`residual_forces`](@ref)
* [`tangent_matrix`](@ref)
* [`strain`](@ref)
* [`stress`](@ref)
* [`structure`](@ref)
"""
abstract type AbstractStructuralState end

Base.copy(x::ST) where {ST<:AbstractStructuralState} = ST([getfield(x, k) for k ∈ fieldnames(ST)]...)

"Returns the `Assembler` used in the `AbstractStructuralState` `st`."
assembler(st::AbstractStructuralState) = st.assembler

"Returns current displacements vector at the current `AbstractStructuralState` `st`."
displacements(st::AbstractStructuralState) = st.Uᵏ

"Returns current displacements increment vector at the current `AbstractStructuralState` `st`."
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Returns current `ResidualsIterationStep` object form an `AbstractStructuralState` `st`."
function iteration_residuals(st::AbstractStructuralState) end

"Returns the current internal forces vector in the `AbstractStructuralState` `st`."
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Returns external forces vector in the `AbstractStructuralState` `st`."
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Returns stresses for each `Element` in the `AbstractStructuralState` `st`."
stress(st::AbstractStructuralState) = st.σᵏ

"Returns strains for each `Element` in the `AbstractStructuralState` `st`."
strain(st::AbstractStructuralState) = st.ϵᵏ

"Returns the structure of the `AbstractStructuralState` `st`."
structure(st::AbstractStructuralState) = st.s

"Returns residual forces vector in the `AbstractStructuralState` `st`."
function residual_forces(st::AbstractStructuralState) end

"Returns system tangent matrix in the `AbstractStructuralState` `st`."
function tangent_matrix(st::AbstractStructuralState) end

"Returns current residuals (forces and displacements) in the `AbstractStructuralState` `st`."
function _residuals(st::AbstractStructuralState)
    RHS_norm = residual_forces(st) |> norm
    Fext_norm = external_forces(st) |> norm
    residual_forces = RHS_norm / Fext_norm

    ΔU_norm = Δ_displacements(st) |> norm
    U_norm = displacements(st) |> norm
    residual_ΔU = ΔU_norm / U_norm
    return residual_forces, residual_ΔU
end


""" Abstract supertype for all structural analysis.

An `AbstractStructuralAnalysis` object facilitates the process of defining an structural analysis
to be solved.

**Common methods:**

* [`structure`](@ref)
* [`free_dofs`](@ref)
* [`structural_state`](@ref)

* [`initial_time`](@ref)
* [`current_time`](@ref)
* [`final_time`](@ref)

* [`current_state`](@ref)
* [`current_iteration`](@ref)
* [`_next!`](@ref)
* [`is_done`](@ref)

"""
abstract type AbstractStructuralAnalysis end

"Returns analyzed structure in the `AbstractStructuralAnalysis` `a`."
structure(a::AbstractStructuralAnalysis) = a.s

"Returns the free dofs of the structure in the `AbstractStructuralAnalysis` `a`."
free_dofs(a::AbstractStructuralAnalysis) = free_dofs(structure(a))

"Returns the initial time of `AbstractStructuralAnalysis` `a`."
initial_time(a::AbstractStructuralAnalysis) = a.t₁

"Returns the current time of `AbstractStructuralAnalysis` `a`."
current_time(a::AbstractStructuralAnalysis) = a.t

"Returns the final time of `AbstractStructuralAnalysis` `a`."
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Increments the time step given the `AbstractStructuralAnalysis` `a` and the `AbstractStructuralSolver` `alg`."
_next!(a::AbstractStructuralAnalysis, alg::AbstractSolver) = a.t += time_step(a)

"Returns `true` if the `AbstractStructuralAnalysis` `a` is completed."
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Returns the current state of the `AbstractStructuralAnalysis` `a`."
current_state(a::AbstractStructuralAnalysis) = a.state

"Returns the current displacements iteration state of the `AbstractStructuralAnalysis` `a`."
current_iteration(a::AbstractStructuralAnalysis) = a.state.iter_state

# ================
# Common methods
# ================

"Applies a fixed displacement boundary condition to the structural analysis `sa` at the current analysis time `t`"
function _apply!(sa::AbstractStructuralAnalysis, lbc::AbstractLoadBoundaryCondition)

    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))

    # Extract dofs to apply the bc
    lbc_dofs_symbols = dofs(lbc)

    # Extract nodes and elements 
    entities = bcs[lbc]
    dofs_lbc = Dof[]

    for dof_symbol in lbc_dofs_symbols
        dofs_lbc_symbol = row_vector(getindex.(dofs.(entities), dof_symbol))
        push!(dofs_lbc, dofs_lbc_symbol...)
    end

    # Repeat the bc values vector to fill a vector of dofs
    dofs_values = lbc(t)
    repeat_mod = Int(length(dofs_lbc) / length(dofs_values))

    external_forces(current_state(sa))[dofs_lbc] = repeat(dofs_values, outer=repeat_mod)

end

"Applies a vector of load boundary conditions to the structure `s` "
_apply!(sa::AbstractStructuralAnalysis, l_bcs::Vector{<:AbstractLoadBoundaryCondition}) = [_apply!(sa, lbc) for lbc in l_bcs]

# ================
# Static analysis
# ================
include("./StaticAnalyses.jl")
@reexport using .StaticAnalyses

# ================
# Dynamic analysis
# ================

# ================
# Modal analysis
# ================

end #module 