"""
Module defining linear static analysis.
This file inherits from src/StructuralAnalyses/StaticAnalyses.jl and extends the methods for
geometric linear static analysis.
"""
module LinearStaticAnalyses

using IterativeSolvers
using Reexport

using ..Utils
using ..StructuralBoundaryConditions
using ..Structures
using ..Assemblers
using ..StaticStates
using ..StructuralAnalyses
using ..StaticAnalyses
using ..StructuralSolvers
using ..Solvers
using ..Solutions

@reexport import ..StructuralSolvers: _solve!, step!

# Since linear analysis do not iterate, the iteration state is:
const LinearResidualsIterationStep = ResidualsIterationStep(nothing, nothing, nothing, nothing, 1,
                                                            ΔU_and_ResidualForce_Criteria())

export LinearStaticAnalysis

"""
A linear analysis is a collection of parameters for defining the static analysis of the structure.
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
As this analysis is linear the stiffness of the structure remains constant at each displacements iteration step.
"""
mutable struct LinearStaticAnalysis{S<:AbstractStructure,R<:Real,LFV<:Vector{R}} <:
               AbstractStaticAnalysis
    "Structure to be analyzed."
    const s::S
    "Structural state."
    const state::FullStaticState
    "Load factors vector of the analysis."
    const λᵥ::LFV
    "Current load factor step."
    current_step::Int64
end

"Constructor for linear analysis with load factors, optional initial step and initial state."
function LinearStaticAnalysis(s::S, λᵥ::LFV;
                              initial_state::FullStaticState=FullStaticState(s,
                                                                             LinearResidualsIterationStep),
                              initial_step::Int=1) where {S<:AbstractStructure,
                                                          LFV<:Vector{<:Real}}
    !(1 ≤ initial_step ≤ length(λᵥ)) &&
        throw(ArgumentError("initial_step must be in [1, $(length(λᵥ))] but is: $initial_step."))
    LinearStaticAnalysis(s, initial_state, λᵥ, initial_step)
end

"Constructor for linear analysis given a final time (or load factor) and the number of steps."
function LinearStaticAnalysis(s::AbstractStructure, final_time::Real=1.0; NSTEPS=10,
                              initial_state::FullStaticState=FullStaticState(s,
                                                                             LinearResidualsIterationStep),
                              initial_step::Int=1)
    t₀ = final_time / NSTEPS
    λᵥ = collect(LinRange(t₀, final_time, NSTEPS))
    LinearStaticAnalysis(s, λᵥ; initial_state, initial_step)
end

function Base.show(io::IO, sa::LinearStaticAnalysis)
    println("LinearStaticAnalysis for:")
    println("• Current load factor $(sa.current_step))/$(length(load_factors(sa))).")
    show(io, sa.s)
    show(io, sa.state)
end

"Solves a linear analysis problem mutating the state."
function _solve!(sa::LinearStaticAnalysis, ::DummySolver)
    s = structure(sa)

    # Initialize solution
    solution = Solution(sa, NewtonRaphson())

    # Load factors iteration.
    while !is_done(sa)

        # Compute external force
        external_forces(current_state(sa)) .= 0
        @debugtime "Assemble external forces" apply!(sa, load_bcs(boundary_conditions(s)))

        # Assemble K
        @debugtime "Assemble internal forces" assemble!(s, sa)

        # Increment structure displacements U = ΔU
        @debugtime "Step" step!(sa)

        # Recompute σ and ε for the assembler
        @debugtime "Assemble internal forces" assemble!(s, sa)

        # Save current state
        push!(solution, current_state(sa))

        # Increments the time or load factor step
        next!(sa)
    end

    solution
end

"Computes ΔU for solving the linear analysis."
function step!(sa::LinearStaticAnalysis)
    # Extract state info
    state = current_state(sa)
    free_dofs_idx = free_dofs(state)

    # Compute Δu
    fₑₓₜ_red = view(external_forces(state), free_dofs_idx)
    K = tangent_matrix(state)[free_dofs_idx, free_dofs_idx]
    ΔU = Δ_displacements(state)
    # NOTE: May not warn if K is singular, see ONSAS#392.
    cg!(ΔU, K, fₑₓₜ_red)

    # Update displacements into the state.
    displacements(state)[free_dofs_idx] .= ΔU
end

end # module
