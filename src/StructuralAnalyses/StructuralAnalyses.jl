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
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
import ..Elements: displacements, external_forces, internal_forces

export AbstractStructuralState, assembler, displacements,
    internal_forces, external_forces, strains, stresses, residual_forces, systemΔu_matrix
export AbstractStructuralAnalysis, structure, initial_time, current_time, final_time,
    _next!, is_done, current_state

""" Abstract supertype for all structural analysis.

An `AbstractStructuralAnalysis` object facilitates the process of defining an structural analysis
to be solved.

**Common methods:**

* [`structure`](@ref)
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
function _next!(a::AbstractStructuralAnalysis, alg::AbstractStructuralSolver) end

"Returns `true` if the `AbstractStructuralAnalysis` `a` is completed."
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Returns the current state of the `AbstractStructuralAnalysis` `a`."
current_state(a::AbstractStructuralAnalysis) = a.state

"Returns the current displacements iteration state of the `AbstractStructuralAnalysis` `a`."
current_iteration(a::AbstractStructuralAnalysis) = a.state.iter_state

# ======================
# Structural state
# ======================

""" Abstract supertype to define a new structural state.
**Common methods:**
* [`assembler`](@ref)
* [`displacements`](@ref)
* [`Δ_displacements`](@ref)
* [`external_forces`](@ref)
* [`internal_forces`](@ref)
* [`residual_forces`](@ref)
* [`tangent_matrix`](@ref)
* [`_tolerancesΔu`](@ref)
* [`strains`](@ref)
* [`stresses`](@ref)
"""
abstract type AbstractStructuralState end

"Returns the `Assembler` used in the `AbstractStructuralState` `st`."
assembler(st::AbstractStructuralState) = st.assembler

"Returns current displacements vector at the current `AbstractStructuralState` `st`."
displacements(st::AbstractStructuralState) = st.Uᵏ

"Returns current displacements increment vector at the current `AbstractStructuralState` `st`."
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Returns the current internal forces vector in the `AbstractStructuralState` `st`."
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Returns external forces vector in the `AbstractStructuralState` `st`."
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Returns stresses for each `Element` in the `AbstractStructuralState` `st`."
stresses(st::AbstractStructuralState) = st.σᵏ

"Returns strains for each `Element` in the `AbstractStructuralState` `st`."
strain(st::AbstractStructuralState) = st.ϵᵏ

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


# ================
# Common methods
# ================
include("./../core/bcs_processing.jl")
include("./../core/assembler.jl")

# ================
# Static analysis
# ================
include("./../analyses/StaticAnalysis.jl")

# ================
# Dynamic analysis
# ================

# ================
# Modal analysis
# ================

end #module 