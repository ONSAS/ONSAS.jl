using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC

using IterativeSolvers: cg

using ..Elements: local_dofs
using ..StructuralModel: AbstractStructure, num_dofs, num_elements
using ..StructuralAnalyses
using ..StructuralSolvers

import ..Utils: _unwrap
import ..StructuralAnalyses: _assemble!, initial_time, current_time, final_time, _next!, residual_forces, tangent_matrix
import ..StructuralSolvers: _init, _solve, _update!

export StaticState
export StaticAnalysis, load_factors, current_load_factor

#==============#
# Static State
#==============#
"""
An `StaticState` object facilitates the process of storing the relevant static variables of the structure. 
### Fields:
- `ΔUᵏ`   -- stores displacements vector increment.
- `Uᵏ`    -- stores displacements vector.
- `Fₑₓₜᵏ` -- stores external forces vector.
- `Fᵢₙₜᵏ` -- stores internal forces vector.
- `assembler object`   -- assembler handler object 
- `iter_state`   -- current Δu iteration state 
"""
mutable struct StaticState <: AbstractStructuralState
    ΔUᵏ::AbstractVector
    Uᵏ::AbstractVector
    Fₑₓₜᵏ::AbstractVector
    Fᵢₙₜᵏ::AbstractVector
    Kₛᵏ::AbstractMatrix
    ϵᵏ::AbstractVector
    σᵏ::AbstractVector
    assembler::Assembler
    iter_state::IterationStep
end

"Returns a default static case for a given mesh."
function StaticState(s::AbstractStructure)
    n_dofs = num_dofs(s)
    n_elements = num_elements(s)
    Uᵏ = @MVector zeros(n_dofs)
    ΔUᵏ = similar(Uᵏ)
    Fₑₓₜᵏ = @MVector zeros(n_dofs)
    Fᵢₙₜᵏ = similar(Uᵏ)
    Kₛᵏ = SparseMatrixCSC(zeros(n_dofs, n_dofs))
    ϵᵏ = Vector{}(undef, n_elements)
    σᵏ = Vector{}(undef, n_elements)
    assemblerᵏ = Assembler(s)
    StaticState(ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assemblerᵏ, IterationStep())
end

residual_forces(sc::StaticState) = sc.Fₑₓₜᵏ - sc.Fᵢₙₜᵏ
tangent_matrix(sc::StaticState) = sc.Kₛᵏ

_unwrap(sc::StaticState) = (sc.ΔUᵏ, sc.Uᵏ, sc.Fₑₓₜᵏ, sc.Fᵢₙₜᵏ, sc.Kₛᵏ, sc.ϵᵏ, sc.σᵏ, sc.assemblerᵏ, sc.iter_state)

#================#
# Static Analysis
#================#
""" StaticAnalysis struct.
A `StaticAnalysis` is a collection of parameters for defining the static analysis of the structure. 
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
### Fields:
- `s`             -- Stores the structure to be analyzed.
- `state`         -- Stores the structural state.
- `λᵥ`            -- Stores the load factors vector of the analysis
- `current_step`  -- Stores the current load factor step
"""
mutable struct StaticAnalysis <: AbstractStructuralAnalysis
    s::AbstractStructure
    state::StaticState
    λᵥ::Vector{<:Real}
    current_step::Int
end

function StaticAnalysis(s::AbstractStructure, λᵥ::Vector{<:Real}; initial_step::Int=0)
    StaticAnalysis(s, StaticState(s), λᵥ, initial_step)
end

function StaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, initial_step::Int=1, init_state::StaticState=StaticState(s))
    t₀ = t₁ / NSTEPS
    λᵥ = LinRange(t₀, t₁, NSTEPS) |> collect
    StaticAnalysis(s, λᵥ, initial_step=initial_step)
end

initial_time(sa::StaticAnalysis) = first(load_factors(sa))

current_time(sa::StaticAnalysis) = load_factors(sa)[sa.current_step]

final_time(sa::StaticAnalysis) = last(load_factors(sa))

"Returns load factors vector"
load_factors(sa::StaticAnalysis) = sa.λᵥ

"Returns the current load factor"
current_load_factor(sa::StaticAnalysis) = current_time(sa)

function _next!(sa::StaticAnalysis, alg::AbstractSolver)
    next_step = sa.current_step + 1
    if next_step > length(load_factors(sa))
        throw(ArgumentError("Analysis is done."))
    else
        sa.current_step = next_step
    end
end

#================#
# Solve
#================#
"Returns the initialized analysis. "
function _init(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    _apply_fixed_bc!(s, sa)

    _update_load_bcs!(s, sa)

    return sa
end

"Internal function to solve different analysis problem"
function _solve(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)


    s = structure(sa)

    # load factors iteration 
    while !is_done(sa)

        _reset!(current_iteration(sa))

        _update_load_bcs!(s, sa)

        @show external_forces(current_state(sa))

        while !_has_converged!(current_iteration(sa), tolerances(alg))

            Main.@infiltrate

            # Computes system residual forces tangent system matrix    
            _assemble_system!(s, sa, alg)

            @show residual_forces(current_state(sa))[s.free_dofs]
            @show tangent_matrix(current_state(sa))[index.(s.free_dofs), index.(s.free_dofs)]

            # Increment U
            @show _step!(sa, alg)

            Main.@infiltrate

        end

        _next!(sa, alg)


    end


    return 1
end

"Rests the assembler state"
function _reset_assembler!(state::StaticState)
    _reset!(state.assembler)
    internal_forces(state) .= 0.0
end



"Computes system residual forces and tangent system matrix for the analysis"
function _assemble_system!(s::AbstractStructure, sa::StaticAnalysis, ::NewtonRaphson)

    s = structure(sa)
    state = current_state(sa)

    # Reset assembler
    _reset_assembler!(state)

    for (mat, mat_elements) in pairs(materials(s))
        for e in mat_elements

            # Global dofs of the element (dofs where K must be added)

            fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e = internal_forces(mat, e, state.Uᵏ)

            # Assembles the element internal magnitudes 
            _assemble!(state, fᵢₙₜ_e, Kᵢₙₜ_e, σ_e, ϵ_e, e)

        end

    end

    # Insert values in the assembler objet into the tangent stiffness matrix
    state.Kₛᵏ = end_assemble(assembler(state))


end

"Assembles the internal force vector and the tangent matrix"
function _assemble!(
    state::StaticState, fᵢₙₜ_e::AbstractVector,
    Kᵢₙₜ_e::AbstractMatrix, σ_e::Real, ϵ_e::Real, e::AbstractElement,
)

    #Returns the local dofs in the global mesh
    ldofs = local_dofs(e)

    # Assembles the element stiffness matrix
    state.Fᵢₙₜᵏ[ldofs] += fᵢₙₜ_e
    push!(state.σᵏ, σ_e)
    push!(state.ϵᵏ, ϵ_e)

    _assemble!(assembler(state), ldofs, Kᵢₙₜ_e)

end


function _step!(sa::StaticAnalysis, alg::NewtonRaphson)

    # Extract state info
    c_state = current_state(sa)
    f_dofs = index.(free_dofs(structure(sa)))

    # Compute Δu
    r = view(residual_forces(c_state), f_dofs)
    K = view(tangent_matrix(c_state), f_dofs, f_dofs)
    @show Δu = cg(K, r)

    # Update sate
    Δ_displacements(c_state)[f_dofs] = Δu
    @show displacements(c_state)[f_dofs] += Δ_displacements(c_state)[f_dofs]

    # Update iteration 
    i_step = _update!(current_iteration(sa),
        Δ_displacements(c_state), displacements(c_state),
        residual_forces(c_state), external_forces(c_state),
        tolerances(alg)
    )

end