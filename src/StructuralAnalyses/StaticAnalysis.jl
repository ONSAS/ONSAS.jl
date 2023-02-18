using StaticArrays: @MVector
using SparseArrays: SparseMatrixCSC

using IterativeSolvers: cg

using ..Elements: local_dofs
using ..StructuralModel: AbstractStructure, num_dofs, num_elements, boundary_conditions, fixed_dof_bcs, load_bcs
using ..StructuralAnalyses
using ..StructuralSolvers

import ..Utils: _unwrap
import ..StructuralAnalyses: _apply!, _assemble!, initial_time, current_time, final_time, _next!, residual_forces, tangent_matrix
import ..StructuralSolvers: _solve, _update!

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
- `assembler`   -- assembler handler object 
- `iter_state`   -- current Δu iteration state 
"""
mutable struct StaticState{DU<:AbstractVector,U<:AbstractVector,
    F<:AbstractVector,K<:AbstractMatrix,E<:AbstractVector,S<:AbstractVector} <: AbstractStructuralState
    ΔUᵏ::DU
    Uᵏ::U
    Fₑₓₜᵏ::F
    Fᵢₙₜᵏ::F
    Kₛᵏ::K
    ϵᵏ::E
    σᵏ::S
    assembler::Assembler
    iter_state::IterationStep
end

# Puede no ser mutable 
# Tengo un strain y stress cambiados en algún lado
"Returns a default static case for a given mesh."
function StaticState(s::AbstractStructure)
    n_dofs = num_dofs(s)
    n_fdofs = num_free_dofs(s)
    n_elements = num_elements(s)
    Uᵏ = @MVector zeros(n_dofs)
    ΔUᵏ = @MVector zeros(n_fdofs)
    Fₑₓₜᵏ = @MVector zeros(n_dofs)
    Fᵢₙₜᵏ = similar(Fₑₓₜᵏ)
    Kₛᵏ = SparseMatrixCSC(zeros(1, 1))
    ϵᵏ = Vector{}(undef, n_elements)
    σᵏ = Vector{}(undef, n_elements)
    assemblerᵏ = Assembler(s)
    StaticState(ΔUᵏ, Uᵏ, Fₑₓₜᵏ, Fᵢₙₜᵏ, Kₛᵏ, ϵᵏ, σᵏ, assemblerᵏ, IterationStep())
end

residual_forces(sc::StaticState, free_dofs::Vector{Dof}) = external_forces(sc)[free_dofs] - internal_forces(sc)[free_dofs]
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
mutable struct StaticAnalysis{S<:AbstractStructure} <: AbstractStructuralAnalysis
    s::S
    state::StaticState
    λᵥ::Vector{<:Real}
    current_step::Int
    function StaticAnalysis(s::S, λᵥ::Vector{<:Real}; initial_step::Int=0) where {S<:AbstractStructure}
        new{S}(s, StaticState(s), λᵥ, initial_step)
    end
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

"Internal function to solve different analysis problem"
function _solve(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)

    states = []

    # load factors iteration 
    while !is_done(sa)

        _reset!(current_iteration(sa)) # Arrancarlo en Infinito

        _apply!(sa, load_bcs(s)) # Compute Fext

        @debug external_forces(current_state(sa))


        while !_has_converged!(current_iteration(sa), tolerances(alg))

            # Computes system residual forces tangent system matrix    
            _assemble_system!(s, sa, alg)

            @debug internal_forces(current_state(sa), s.free_dofs)
            @debug residual_forces(current_state(sa), s.free_dofs)
            @debug tangent_matrix(current_state(sa))[index.(s.free_dofs), index.(s.free_dofs)]

            # Increment U 
            @debug _step!(sa, alg)

        end

        push!(states, current_state(sa))

        _next!(sa, alg)

    end


    return states
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
    f_dofs = free_dofs(structure(sa))
    f_dofs_indexes = index.(f_dofs)

    # Compute Δu
    r = residual_forces(c_state, f_dofs)
    K = view(tangent_matrix(c_state), f_dofs_indexes, f_dofs_indexes)
    Δ_displacements(c_state) = cg(K, r)

    # Update sate
    displacements(c_state)[f_dofs] += Δ_displacements(c_state)

    # Update iteration 
    i_step = _update!(current_iteration(sa),
        Δ_displacements(c_state), displacements(c_state),
        residual_forces(c_state, f_dofs), external_forces(c_state),
        tolerances(alg)
    )

end