"""
Module defining structural analyses that can be solved.
Each structural analysis consists of a data type with an structure to be analyzed,
analysis parameters and an specific state which describes the current state of structure
in the analysis.
"""
module StructuralAnalyses

using LinearAlgebra, Reexport

using ..Materials
using ..Entities
using ..BoundaryConditions
using ..Meshes
using ..Structures
using ..StructuralSolvers
using ..Utils

@reexport import ..Utils: apply!
@reexport import ..Entities: internal_forces, inertial_forces, strain, stress
@reexport import ..Structures: free_dofs
@reexport import ..StructuralSolvers: _update!
@reexport import ..Assemblers: assemble!, _end_assemble!
@reexport import ..StructuralSolvers: reset!
@reexport import ..Solutions: displacements, external_forces, iteration_residuals

export AbstractStructuralState, Δ_displacements, tangent_matrix, residual_forces!, tangent_matrix,
       structure, assembler, residual_forces_norms, residual_displacements_norms,
       AbstractStructuralAnalysis, initial_time, current_time, final_time, _next!, is_done,
       current_state, current_iteration

""" Abstract supertype to define a new structural state.
**Abstract Methods**
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
* [`assemble!`](@ref)
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

"Return free `Dof`s of the structure in the `AbstractStructuralState` `st`."
free_dofs(st::AbstractStructuralState) = st.free_dofs

# Assemble
"Assembles the element `e` internal forces `fᵢₙₜ_e` into the `AbstractState` `st`"
function assemble!(st::AbstractStructuralState, fᵢₙₜ_e::AbstractVector, e::AbstractElement)
    view(internal_forces(st), local_dofs(e)) .+= fᵢₙₜ_e
end

"Assembles the element `e` stiffness matrix matrix `K_e` into the `AbstractState` `st`"
function assemble!(st::AbstractStructuralState, kₛ_e::AbstractMatrix, e::AbstractElement)
    assemble!(assembler(st), local_dofs(e), kₛ_e)
end

"Assembles the element `e` stress σₑ and strain ϵₑ into the `AbstractState` `st`"
function assemble!(st::AbstractStructuralState, σₑ::E, ϵₑ::E,
                   e::AbstractElement) where {E<:Union{Real,AbstractMatrix}}
    stress(st)[e] .= σₑ
    strain(st)[e] .= ϵₑ
end

"Fill the system tangent matrix in the `AbstractStructuralState` `st` once the `Assembler` object is built."
_end_assemble!(st::AbstractStructuralState) = _end_assemble!(tangent_matrix(st), assembler(st))

"Return system tangent matrix in the `AbstractStructuralState` `st`."
function tangent_matrix(st::AbstractStructuralState, alg::AbstractSolver) end

"Return relative residual forces for the current `AbstractStructuralState` `st`."
function residual_forces_norms(st::AbstractStructuralState)
    rᵏ_norm = norm(residual_forces!(st))
    fₑₓₜ_norm = norm(external_forces(st))
    rᵏ_norm, rᵏ_norm / fₑₓₜ_norm
end

"Return relative residual displacements for the current `AbstractStructuralState` `st`."
function residual_displacements_norms(st::AbstractStructuralState)
    ΔU_norm = norm(Δ_displacements(st))
    U_norm = norm(displacements(st))
    ΔU_norm, ΔU_norm / U_norm
end

"Updates the `AbstractStructuralState` `st` during the displacements iteration."
function _update!(st::AbstractStructuralState, args...; kwargs...) end

"Resets  the `AbstractStructuralState` assembled magnitudes before starting a new assembly."
function reset!(st::AbstractStructuralState, args...; kwargs...) end

""" Abstract supertype for all structural analysis.

An structural analysis object facilitates the process of defining an structural analysis
to be solved.

**Abstract Methods**

* [`structure`](@ref)
* [`free_dofs`](@ref)

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

"Return analyzed structure in the structural analysis."
structure(a::AbstractStructuralAnalysis) = a.s

"Return the initial time of structural analysis."
initial_time(a::AbstractStructuralAnalysis) = a.t₁

"Return the current time of structural analysis."
current_time(a::AbstractStructuralAnalysis) = a.t

"Return the final time of structural analysis."
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Increment the time step given of a structural analysis. Dispatch is done for different
solvers."
_next!(a::AbstractStructuralAnalysis, solver::AbstractSolver) = a.t += time_step(a)

"Return true if the structural analysis is completed."
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Return the current state of the structural analysis."
current_state(a::AbstractStructuralAnalysis) = a.state

"Return the current displacements iteration state of the structural analysis."
current_iteration(a::AbstractStructuralAnalysis) = iteration_residuals(a.state)

"Rests the structural analysis (sets the current time to the initial time)."
function reset!(a::AbstractStructuralState) end

# ================
# Common methods
# ================

"Apply a boundary condition to the structural analysis at the current analysis time."
function apply!(sa::AbstractStructuralAnalysis, lbc::AbstractNeumannBoundaryCondition)
    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))
    dofs_lbc, dofs_values = apply(bcs, lbc, t)
    external_forces(current_state(sa))[dofs_lbc] = dofs_values
end

"Apply a vector of load boundary conditions to the structure."
function apply!(sa::AbstractStructuralAnalysis, l_bcs::Vector{<:AbstractNeumannBoundaryCondition})
    [apply!(sa, lbc) for lbc in l_bcs]
end

# ================
# Dynamic analysis
# ================

# ================
# Modal analysis
# ================

end #module
