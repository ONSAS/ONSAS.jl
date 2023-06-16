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

@reexport import ..StructuralSolvers: _solve!, _step!

export LinearStaticAnalysis

"""
A linear analysis is a collection of parameters for defining the static analysis of the structure.
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
As this analysis is linear the stiffness of the structure remains constant at each displacements iteration step.
"""
struct LinearStaticAnalysis{S<:AbstractStructure,LFV<:AbstractVector{<:Real}} <:
       AbstractStaticAnalysis
    "Structure to be analyzed."
    s::S
    "Structural state."
    state::StaticState
    "Load factors vector of the analysis."
    λᵥ::LFV
    "Current load factor step."
    current_step::ScalarWrapper{Int}
    function LinearStaticAnalysis(s::S, λᵥ::LFV;
                                  initial_step::Int=1) where {S<:AbstractStructure,
                                                              LFV<:AbstractVector{<:Real}}
        # Since linear analysis is not iterating
        iter_state = ResidualsIterationStep(nothing, nothing, nothing, nothing, 0,
                                            ΔU_and_ResidualForce_Criteria())
        new{S,LFV}(s, StaticState(s, iter_state), λᵥ, ScalarWrapper(initial_step))
    end
end

"Constructor for linear analysis given a final time (or load factor) and the number of steps."
function LinearStaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, initial_step::Int=1)
    t₀ = t₁ / NSTEPS
    λᵥ = collect(LinRange(t₀, t₁, NSTEPS))
    return LinearStaticAnalysis(s, λᵥ; initial_step=initial_step)
end

function Base.show(io::IO, sa::LinearStaticAnalysis)
    println("LinearStaticAnalysis for:")
    show(io, sa.s)
    show(io, sa.state)
end

"Solves a linear analysis problem mutating the state."
function _solve!(sa::LinearStaticAnalysis)
    s = structure(sa)

    # Initialize solution
    solution = StatesSolution(sa, NewtonRaphson())

    # Load factors iteration.
    while !is_done(sa)
        # Set displacements to zero
        displacements(current_state(sa)) .= 0.0

        # Compute external force
        apply!(sa, load_bcs(boundary_conditions(s))) # Compute Fext

        # Assemble K
        _assemble!(s, sa)

        # Increment structure displacements U = U + ΔU
        _step!(sa)

        # Recompute σ and ε for the assembler
        _assemble!(s, sa)

        # Save current state
        push!(solution, current_state(sa))

        # Increments the time or load factor step
        _next!(sa)
    end

    solution
end

"Computes ΔU for solving the linear analysis."
function _step!(sa::LinearStaticAnalysis)
    # Extract state info
    state = current_state(sa)
    free_dofs_idx = free_dofs(state)

    # Compute Δu
    fₑₓₜ_red = view(external_forces(state), free_dofs_idx)
    K = tangent_matrix(state)[free_dofs_idx, free_dofs_idx]
    ΔU = Δ_displacements(state)
    cg!(ΔU, K, fₑₓₜ_red)

    # Update displacements into the state.
    state.Uᵏ[free_dofs_idx] .+= ΔU
end

end # module
