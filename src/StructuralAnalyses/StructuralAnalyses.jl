"""
Module defining structural analyses that can be solved. 
Each structural analysis consists of a data type with an structure to be analyzed, analysis parameters and 
an specific state which describes the current state of structure in the analysis.
"""
module StructuralAnalyses

using LinearAlgebra: norm
using Reexport: @reexport

@reexport using ..Materials
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
@reexport using ..StructuralModel
@reexport using ..StructuralSolvers
using ..Utils: row_vector

import ..Elements: internal_forces, inertial_forces, strain, stress
import ..StructuralModel: free_dofs
import ..StructuralSolvers: _assemble!, _update!, _reset!, _end_assemble!

export AbstractStructuralState, _apply!, _assemble!, displacements, Δ_displacements, external_forces, residual_forces,
    tangent_matrix, residual_forces_norms, residual_displacements_norms, iteration_residuals, tangent_matrix, structure,
    assembler, iteration_residuals, residual_forces_norms, residual_displacements_norms

export AbstractStructuralAnalysis, initial_time, current_time, final_time, _next!, is_done, current_state, current_iteration

""" Abstract supertype to define a new structural state.
**Common methods:**
### Accessors:
* [`displacements`](@ref)
* [`Δ_displacements`](@ref)
* [`external_forces`](@ref)
* [`internal_forces`](@ref)
* [`residual_forces`](@ref)
* [`tangent_matrix`](@ref)
* [`strain`](@ref)
* [`stress`](@ref)
* [`structure`](@ref)
* [`free_dofs`](@ref)


### Iteration:
* [`_assemble!`](@ref)
* [`assembler`](@ref)
* [`iteration_residuals`](@ref)
* [`residual_forces_norms`](@ref)
* [`residual_displacements_norms`](@ref)
"""
abstract type AbstractStructuralState end

#Accessors

"Returns the `Assembler` used in the `AbstractStructuralState` `st`."
assembler(st::AbstractStructuralState) = st.assembler

"Returns current `ResidualsIterationStep` object form an `AbstractStructuralState` `st`."
iteration_residuals(st::AbstractStructuralState) = st.iter_state

"Returns current displacements vector at the current `AbstractStructuralState` `st`."
displacements(st::AbstractStructuralState) = st.Uᵏ

"Returns current displacements increment vector at the current `AbstractStructuralState` `st`."
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Returns the current internal forces vector in the `AbstractStructuralState` `st`."
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Returns external forces vector in the `AbstractStructuralState` `st`."
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Returns residual forces vector in the `AbstractStructuralState` `st`."
function residual_forces(st::AbstractStructuralState) end

"Returns stresses for each `Element` in the `AbstractStructuralState` `st`."
stress(st::AbstractStructuralState) = st.σᵏ

"Returns strains for each `Element` in the `AbstractStructuralState` `st`."
strain(st::AbstractStructuralState) = st.ϵᵏ

"Returns the structure of the `AbstractStructuralState` `st`."
structure(st::AbstractStructuralState) = st.s

"Returns free `Dof`s of the structure in the `AbstractStructuralState` `st`."
free_dofs(st::AbstractStructuralState) = free_dofs(structure(st))

# Assemble
"Assembles the element `e` internal forces `fᵢₙₜ_e` into the `AbstractState` `st`"
_assemble!(st::AbstractStructuralState, fᵢₙₜ_e::AbstractVector, e::AbstractElement) =
    view(internal_forces(st), index.(local_dofs(e))) .+= fᵢₙₜ_e

"Assembles the element `e` stiffness matrix matrix `K_e` into the `AbstractState` `st`"
_assemble!(st::AbstractStructuralState, kₛ_e::AbstractMatrix, e::AbstractElement) =
    _assemble!(assembler(st), local_dofs(e), kₛ_e)

"Assembles the element `e` stress σₑ and strain ϵₑ into the `AbstractState` `st`"
function _assemble!(st::AbstractStructuralState, σₑ::E, ϵₑ::E, e::AbstractElement) where {E<:Union{Real,AbstractMatrix}}
    stress(st)[e] = σₑ
    strain(st)[e] = ϵₑ
end

"Fill the system tangent matrix in the `AbstractStructuralState` `st` once the `Assembler` object is built."
_end_assemble!(st::AbstractStructuralState) = _end_assemble!(tangent_matrix(st), assembler(st))

"Returns system tangent matrix in the `AbstractStructuralState` `st`."
function tangent_matrix(st::AbstractStructuralState, alg::AbstractSolver) end

"Returns relative residual forces for the current `AbstractStructuralState` `st`."
function residual_forces_norms(st::AbstractStructuralState)
    rᵏ_norm = residual_forces(st) |> norm
    fₑₓₜ_norm = external_forces(st) |> norm
    return rᵏ_norm, rᵏ_norm / fₑₓₜ_norm
end

"Returns relative residual displacements for the current `AbstractStructuralState` `st`."
function residual_displacements_norms(st::AbstractStructuralState)
    ΔU_norm = Δ_displacements(st) |> norm
    U_norm = displacements(st) |> norm
    return ΔU_norm, ΔU_norm / U_norm
end

"Updates the `AbstractStructuralState` `st` during the displacements iteration."
function _update!(st::AbstractStructuralState, args...; kwargs...) end

"Resets  the `AbstractStructuralState` assembled magnitudes before starting a new assembly."
function _reset!(st::AbstractStructuralState, args...; kwargs...) end

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
current_iteration(a::AbstractStructuralAnalysis) = iteration_residuals(a.state)

# ================
# Common methods
# ================

"Applies an `AbstractLoadBoundaryCondition` lbc into the structural analysis `sa` at the current analysis time `t`"
function _apply!(sa::AbstractStructuralAnalysis, lbc::AbstractLoadBoundaryCondition)

    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))
    dofs_lbc, dofs_values = _apply(bcs, lbc, t)

    external_forces(current_state(sa))[dofs_lbc] = dofs_values

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