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
using ..Utils

@reexport import ..Utils: apply!
@reexport import ..Entities: internal_forces, inertial_forces, strain, stress, elements_cache
@reexport import ..Structures: free_dofs
@reexport import ..Assemblers: assemble!, end_assemble!

export AbstractStructuralState, AbstractStaticState, AbstractDynamicState, Δ_displacements,
       Δ_displacements!, residual_forces!, structure, assembler, residual_forces_norms,
       residual_displacements_norms, AbstractStructuralAnalysis, initial_time, current_time, times,
       final_time, is_done, current_state, current_iteration, displacements, external_forces,
       iteration_residuals, tangent_matrix, internal_cache, elements_cache, velocity, acceleration,
       viscous_forces, mass_matrix, damping_matrix, stiffness_matrix

""" Abstract supertype to define a new structural state.
**Abstract Methods**
### Accessors:
* [`displacements`](@ref)
* [`Δ_displacements`](@ref)
* [`Δ_displacements!`](@ref)
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

"""
States representing static analyses.
"""
abstract type AbstractStaticState <: AbstractStructuralState end

"Return the assembler used in the structural state."
assembler(st::AbstractStructuralState) = st.assembler

"Return current `ResidualsIterationStep` object form an structural state."
iteration_residuals(st::AbstractStructuralState) = st.iter_state

"Return current displacements vector at the current structural state."
displacements(st::AbstractStructuralState) = st.Uᵏ

"Return current displacements increment vector at the current structural state."
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Update and return current displacements increment vector at the current structural state."
function Δ_displacements!(st::AbstractStructuralState, ΔUᵏ⁺¹::AbstractVector)
    st.ΔUᵏ .= ΔUᵏ⁺¹
    st.ΔUᵏ
end

"Return the current internal forces vector in the structural state."
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Return external forces vector in the structural state."
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Return residual forces vector in the structural state."
function residual_forces!(st::AbstractStructuralState) end

"Return system tangent matrix in the structural state given a solver."
tangent_matrix(st::AbstractStructuralState) = st.Kₛᵏ

"Return stresses for each `Element` in the structural state."
stress(st::AbstractStructuralState) = st.σᵏ

"Return strains for each `Element` in the structural state."
strain(st::AbstractStructuralState) = st.ϵᵏ

"Return free `Dof`s of the structure in the structural state."
free_dofs(st::AbstractStructuralState) = st.free_dofs

# Assemble
"Assembles the element `e` internal forces `fᵢₙₜ_e` into the structural state."
function assemble!(st::AbstractStructuralState, fᵢₙₜ_e::AbstractVector, e::AbstractElement)
    view(internal_forces(st), local_dofs(e)) .+= fᵢₙₜ_e
end

"Assembles the element `e` stiffness matrix matrix `K_e` into the structural state."
function assemble!(st::AbstractStructuralState, kₛ_e::AbstractMatrix, e::AbstractElement)
    assemble!(assembler(st), local_dofs(e), kₛ_e)
end

"Assembles the element `e` stress σₑ and strain ϵₑ into the structural state."
function assemble!(st::AbstractStructuralState, σₑ::ST, ϵₑ::ET,
                   e::AbstractElement) where {ST<:Union{Real,AbstractMatrix},
                                              ET<:Union{Real,AbstractMatrix}}
    stress(st)[e] .= σₑ
    strain(st)[e] .= ϵₑ
end

"Fill the system tangent matrix in the structural state once the assembler object is built."
end_assemble!(st::AbstractStructuralState) = end_assemble!(tangent_matrix(st), assembler(st))

"Return relative residual forces for the current structural state."
function residual_forces_norms(st::AbstractStructuralState)
    rᵏ_norm = norm(residual_forces!(st))
    fₑₓₜ_norm = norm(external_forces(st))
    rᵏ_norm, rᵏ_norm / fₑₓₜ_norm
end

"Return relative residual displacements for the current structural state."
function residual_displacements_norms(st::AbstractStructuralState)
    ΔU_norm = norm(Δ_displacements(st))
    U_norm = norm(displacements(st))
    ΔU_norm, ΔU_norm / U_norm
end

"Return the cache associated to the given element."
function internal_cache(state::AbstractStructuralState, e::AbstractElement)
    internal_cache(state, typeof(e))
end
internal_cache(::AbstractStructuralState, ::Type{<:AbstractElement}) = nothing

function elements_cache(s::AbstractStructuralState, e::AbstractElement)
    elements_cache(assembler(s), e)
end

"""
States representing dynamic analyses.
"""
abstract type AbstractDynamicState <: AbstractStructuralState end

"Return current velocity vector at the current structural state."
velocity(st::AbstractDynamicState) = st.Udotᵏ

"Return current acceleration vector at the current structural state."
acceleration(st::AbstractDynamicState) = st.Udotdotᵏ

"Return the current inertial forces vector in the structural state."
inertial_forces(st::AbstractStructuralState) = st.Fᵢₙₑᵏ

"Return the current viscous forces vector in the structural state."
viscous_forces(st::AbstractStructuralState) = st.Fᵥᵢₛᵏ

"Return the current mass matrix in the structural state."
mass_matrix(st::AbstractStructuralState) = st.Mᵏ

"Return the current damping matrix in the structural state."
damping_matrix(st::AbstractStructuralState) = st.Cᵏ

"Return the current stiffness matrix in the structural state."
stiffness_matrix(st::AbstractStructuralState) = st.Kᵏ

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
* [`next!`](@ref)
* [`is_done`](@ref)
"""
abstract type AbstractStructuralAnalysis end

"Return analyzed structure in the structural analysis."
structure(a::AbstractStructuralAnalysis) = a.s

"Return the initial time of structural analysis."
initial_time(a::AbstractStructuralAnalysis) = a.t₁

"Return the current time of structural analysis."
current_time(a::AbstractStructuralAnalysis) = a.t

"Return the time vector of a structural analysis."
times(a::AbstractStructuralAnalysis) = a.t

"Return the final time of structural analysis."
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Return true if the structural analysis is completed."
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Return the current state of the structural analysis."
current_state(a::AbstractStructuralAnalysis) = a.state

"Return the current displacements iteration state of the structural analysis."
current_iteration(a::AbstractStructuralAnalysis) = iteration_residuals(a.state)

# ================
# Common methods
# ================

"Apply a boundary condition to the structural analysis at the current analysis time."
function apply!(sa::AbstractStructuralAnalysis, lbc::AbstractNeumannBoundaryCondition)
    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))
    dofs_lbc, dofs_values = apply(bcs, lbc, t)
    external_forces(current_state(sa))[dofs_lbc] .+= dofs_values
end

"Apply a vector of load boundary conditions to the structure."
function apply!(sa::AbstractStructuralAnalysis, l_bcs::Vector{<:AbstractNeumannBoundaryCondition})
    for lbc in l_bcs
        apply!(sa, lbc)
    end
end

end
