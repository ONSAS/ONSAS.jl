"""
Module defining structural analyses that can be solved. 
"""
module StructuralAnalyses

using Reexport: @reexport

@reexport using ..StructuralModel
@reexport using ..Materials
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
@reexport import ..Utils: displacements, external_forces, internal_forces

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

"Returns analyzed structure"
structure(a::AbstractStructuralAnalysis) = a.s

"Returns initial time of the analysis"
initial_time(a::AbstractStructuralAnalysis) = a.t₁

"Returns current time of the analysis"
current_time(a::AbstractStructuralAnalysis) = a.t

"Returns final time of the analysis"
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Increments the current time of the analysis"
function _next!(a::AbstractStructuralAnalysis) end

"Returns `true` if the analysis is done"
is_done(a::AbstractStructuralAnalysis) = current_time(a) > final_time(a)

"Returns the current structural state"
current_state(a::AbstractStructuralAnalysis) = a.state

"Returns the current iteration state"
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

"Returns an assembler object (this is used to speed up the assembling process"
assembler(st::AbstractStructuralState) = st.assembler

"Returns current displacement vector"
displacements(st::AbstractStructuralState) = st.Uᵏ

"Returns current displacement increments vector"
Δ_displacements(st::AbstractStructuralState) = st.ΔUᵏ

"Returns current internal forces vector"
internal_forces(st::AbstractStructuralState) = st.Fᵢₙₜᵏ

"Returns current external forces vector"
external_forces(st::AbstractStructuralState) = st.Fₑₓₜᵏ

"Returns current stresses"
stresses(st::AbstractStructuralState) = st.σᵏ

"Returns current strains"
strain(st::AbstractStructuralState) = st.ϵᵏ

"Returns current external forces vector"
function residual_forces(st::AbstractStructuralState) end

"Returns current tangent matrix"
function tangent_matrix(st::AbstractStructuralState) end

"Returns forces and displacements in u tolerances"
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