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

import ..Elements: internal_forces, inertial_forces, strain, stress
import ..StructuralModel: free_dofs
import ..StructuralSolvers: _assemble!, _update!, _end_assemble!
@reexport import ..StructuralSolvers: displacements, external_forces, iteration_residuals

export AbstractStructuralState, _apply!, _assemble!, Δ_displacements, tangent_matrix,
       residual_forces!,
       tangent_matrix, structure, assembler, residual_forces_norms, residual_displacements_norms

export AbstractStructuralAnalysis, initial_time, current_time, final_time, _next!, is_done,
       current_state, current_iteration, reset!

""" Abstract supertype to define a new structural state.
**Common methods:**
### Accessors:
* [`displacements`](@ref)
* [`Δ_displacements`](@ref)
* [`external_forces`](@ref)
* [`internal_forces`](@ref)
* [`residual_forces!`](@ref)
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

"Return the `Assembler` used in the `AbstractStructuralState` `st`."
assembler(st::AbstractStructuralState) = st.assembler

"Return current `ResidualsIterationStep` object form an `AbstractStructuralState` `st`."
iteration_residuals(st::AbstractStructuralState) = st.iter_state

"Return current displacements vector at the current `AbstractStructuralState` `st`."
displacements(st::AbstractStructuralState) = st.Uᵏ

"Return current displacements increment vector at the current `AbstractStructuralState` `st`."
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Return the current internal forces vector in the `AbstractStructuralState` `st`."
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Return external forces vector in the `AbstractStructuralState` `st`."
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Return residual forces vector in the `AbstractStructuralState` `st`."
function residual_forces!(st::AbstractStructuralState) end

"Return stresses for each `Element` in the `AbstractStructuralState` `st`."
stress(st::AbstractStructuralState) = st.σᵏ

"Return strains for each `Element` in the `AbstractStructuralState` `st`."
strain(st::AbstractStructuralState) = st.ϵᵏ

"Return the structure of the `AbstractStructuralState` `st`."
structure(st::AbstractStructuralState) = st.s

"Return free `Dof`s of the structure in the `AbstractStructuralState` `st`."
free_dofs(st::AbstractStructuralState) = free_dofs(structure(st))

# Assemble
"Assembles the element `e` internal forces `fᵢₙₜ_e` into the `AbstractState` `st`"
function _assemble!(st::AbstractStructuralState, fᵢₙₜ_e::AbstractVector, e::AbstractElement)
    return view(internal_forces(st), local_dofs(e)) .+= fᵢₙₜ_e
end

"Assembles the element `e` stiffness matrix matrix `K_e` into the `AbstractState` `st`"
function _assemble!(st::AbstractStructuralState, kₛ_e::AbstractMatrix, e::AbstractElement)
    return _assemble!(assembler(st), local_dofs(e), kₛ_e)
end

"Assembles the element `e` stress σₑ and strain ϵₑ into the `AbstractState` `st`"
function _assemble!(st::AbstractStructuralState, σₑ::E, ϵₑ::E,
                    e::AbstractElement) where {E<:Union{Real,AbstractMatrix}}
    stress(st)[e] .= σₑ
    return strain(st)[e] .= ϵₑ
end

"Fill the system tangent matrix in the `AbstractStructuralState` `st` once the `Assembler` object is built."
_end_assemble!(st::AbstractStructuralState) = _end_assemble!(tangent_matrix(st), assembler(st))

"Return system tangent matrix in the `AbstractStructuralState` `st`."
function tangent_matrix(st::AbstractStructuralState, alg::AbstractSolver) end

"Return relative residual forces for the current `AbstractStructuralState` `st`."
function residual_forces_norms(st::AbstractStructuralState)
    rᵏ_norm = norm(residual_forces!(st))
    fₑₓₜ_norm = norm(external_forces(st))
    return rᵏ_norm, rᵏ_norm / fₑₓₜ_norm
end

"Return relative residual displacements for the current `AbstractStructuralState` `st`."
function residual_displacements_norms(st::AbstractStructuralState)
    ΔU_norm = norm(Δ_displacements(st))
    U_norm = norm(displacements(st))
    return ΔU_norm, ΔU_norm / U_norm
end

"Updates the `AbstractStructuralState` `st` during the displacements iteration."
function _update!(st::AbstractStructuralState, args...; kwargs...) end

"Resets  the `AbstractStructuralState` assembled magnitudes before starting a new assembly."
function reset!(st::AbstractStructuralState, args...; kwargs...) end

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
* [`reset!`](@ref)

"""
abstract type AbstractStructuralAnalysis end

"Return analyzed structure in the `AbstractStructuralAnalysis` `a`."
structure(a::AbstractStructuralAnalysis) = a.s

"Return the initial time of `AbstractStructuralAnalysis` `a`."
initial_time(a::AbstractStructuralAnalysis) = a.t₁

"Return the current time of `AbstractStructuralAnalysis` `a`."
current_time(a::AbstractStructuralAnalysis) = a.t

"Return the final time of `AbstractStructuralAnalysis` `a`."
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Increments the time step given the `AbstractStructuralAnalysis` `a` and the `AbstractStructuralSolver` `alg`."
_next!(a::AbstractStructuralAnalysis, alg::AbstractSolver) = a.t += time_step(a)

"Return `true` if the `AbstractStructuralAnalysis` `a` is completed."
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Return the current state of the `AbstractStructuralAnalysis` `a`."
current_state(a::AbstractStructuralAnalysis) = a.state

"Return the current displacements iteration state of the `AbstractStructuralAnalysis` `a`."
current_iteration(a::AbstractStructuralAnalysis) = iteration_residuals(a.state)

"Rests the `AbstractStructuralAnalysis` `sa` (sets the current time to the initial time)."
function reset!(a::AbstractStructuralState) end

# ================
# Common methods
# ================

"Apply an `AbstractLoadBoundaryCondition` lbc into the structural analysis `sa` at the current analysis time `t`"
function _apply!(sa::AbstractStructuralAnalysis, lbc::AbstractLoadBoundaryCondition)
    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))
    dofs_lbc, dofs_values = _apply(bcs, lbc, t)
    return external_forces(current_state(sa))[dofs_lbc] = dofs_values
end

"Apply a vector of load boundary conditions to the structure `s` "
function _apply!(sa::AbstractStructuralAnalysis, l_bcs::Vector{<:AbstractLoadBoundaryCondition})
    return [_apply!(sa, lbc) for lbc in l_bcs]
end

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
