module StaticAnalyses

using Dictionaries: dictionary
using Reexport: @reexport

@reexport using ..Materials
@reexport using ...Elements
@reexport using ...Meshes
@reexport using ...StructuralModel: AbstractStructure, load_bcs
@reexport using ..StructuralAnalyses
@reexport using ...StructuralSolvers

import ..StructuralAnalyses: _assemble!, initial_time, current_time, final_time, _next!,
    iteration_residuals, is_done, reset!

export AbstractStaticAnalysis, load_factors, current_load_factor

include("StaticState.jl")

""" Abstract supertype for all structural analysis.

An `AbstractStaticAnalysis` object facilitates the process of defining an static analysis
to be solved. Time variable for static analysis is used to obtain a load factor value.  
Of course this abstract type inherits from `AbstractStructuralAnalysis` type, 
and extends the following methods:

**Common methods:**

* [`initial_time`](@ref)
* [`current_time`](@ref)
* [`final_time`](@ref)
* [`load_factors`](@ref)
* [`current_load_factor`](@ref)
* [`_next!`](@ref)
* [`is_done`](@ref)
* [`reset!`](@ref)
* [`_solve`](@ref)
* [`_step!`](@ref)

**Common fields**
- `λᵥ`           -- stores the load factors vector.
- `current_step` -- stores the current step of the analysis.
- `state`        -- stores the current state of the analysis.


"""
abstract type AbstractStaticAnalysis <: AbstractStructuralAnalysis end

"Returns the initial load factor of an `AbstractStaticAnalysis` `sa`."
initial_time(sa::AbstractStaticAnalysis) = first(load_factors(sa))

"Returns the current load factor of an `AbstractStaticAnalysis` `sa`."
current_time(sa::AbstractStaticAnalysis) = load_factors(sa)[sa.current_step[]]

"Returns the final load factor of an `AbstractStaticAnalysis` `sa`."
final_time(sa::AbstractStaticAnalysis) = last(load_factors(sa))

"Returns `true` if the `AbstractStaticAnalysis` `sa` is completed."
function is_done(sa::AbstractStaticAnalysis)
    is_done_bool = if sa.current_step[] > length(load_factors(sa))
        sa.current_step[] -= 1
        true
    else
        false
    end
end

"Returns the final load factor vector of an `AbstractStaticAnalysis` `sa`."
load_factors(sa::AbstractStaticAnalysis) = sa.λᵥ

"Returns the current load factor of an `AbstractStaticAnalysis` `sa`."
current_load_factor(sa::AbstractStaticAnalysis) = current_time(sa)

"Jumps to the next current load factor defined in the `AbstractStaticAnalysis` `sa`."
_next!(sa::AbstractStaticAnalysis) = sa.current_step[] += 1

"Sets the current load factor of the `AbstractStaticAnalysis` `sa` to the initial load factor.
Also resets! the iteration and `AbstractStructuralState`."
function reset!(sa::AbstractStaticAnalysis)
    sa.current_step[] = 1
    reset!(current_state(sa))
    @info "The current time of analysis have been reset."
    sa
end

"Assembles the Structure `s` (internal forces) during the `StaticAnalysis` `sa`."
function _assemble!(s::AbstractStructure, sa::AbstractStaticAnalysis)

    state = current_state(sa)

    # Reset assembled magnitudes
    _reset_assemble!(state)

    for (mat, mat_elements) in pairs(materials(s))
        for e in mat_elements

            # Global dofs of the element (dofs where K must be added)
            u_e = view(displacements(state), local_dofs(e))
            fᵢₙₜ_e, kₛ_e, σ_e, ϵ_e = internal_forces(mat, e, u_e)

            # Assembles the element internal magnitudes 
            _assemble!(state, fᵢₙₜ_e, e)
            _assemble!(state, kₛ_e, e)
            _assemble!(state, σ_e, ϵ_e, e)

        end

    end

    # Insert values in the assembler objet into the sysyem tangent stiffness matrix
    _end_assemble!(state)

end

"Resets the assembled magnitudes of the `AbstractStructuralState` `state`."
function _reset_assemble!(state::AbstractStructuralState)
    _reset!(assembler(state))
    internal_forces(state) .= 0.0
    tangent_matrix(state)[findall(!iszero, tangent_matrix(state))] .= 0.0
    # FIXME: Zero out stress and strain
    nothing
end

"Pushes the current state `c_state` into the `StatesSolution` `st_sol`."
function Base.push!(st_sol::StatesSolution, c_state::StaticState)

    # Pointers
    s = structure(c_state)
    # Empty assembler since the info is stored in k
    assemblerᵏ = Assembler(s)

    # Deep copies 
    Uᵏ = deepcopy(displacements(c_state))
    ΔUᵏ = deepcopy(Δ_displacements(c_state))
    fₑₓₜᵏ = deepcopy(external_forces(c_state))
    fᵢₙₜᵏ = deepcopy(internal_forces(c_state))
    Kₛᵏ = deepcopy(tangent_matrix(c_state))
    res_forces = deepcopy(c_state.res_forces)
    σᵏ = dictionary([e => deepcopy(σ) for (e, σ) in pairs(stress(c_state))])
    ϵᵏ = dictionary([e => deepcopy(ϵ) for (e, ϵ) in pairs(strain(c_state))])
    iter_state = deepcopy(iteration_residuals(c_state))

    push!(states(st_sol), StaticState(s, ΔUᵏ, Uᵏ, fₑₓₜᵏ, fᵢₙₜᵏ, Kₛᵏ, res_forces, ϵᵏ, σᵏ, assemblerᵏ, iter_state))
end

include("LinearStaticAnalyses.jl")
@reexport using .LinearStaticAnalyses

include("NonLinearStaticAnalyses.jl")
@reexport using .NonLinearStaticAnalyses

end # module
