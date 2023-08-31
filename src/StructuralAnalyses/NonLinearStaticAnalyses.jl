"""
Module defining non-linear static analysis.
This file inherits from src/StructuralAnalyses/StaticAnalyses.jl and extends the methods for
geometric non-linear static analysis.
"""
module NonLinearStaticAnalyses

using LinearAlgebra: norm
using IterativeSolvers: cg!
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

export NonLinearStaticAnalysis

"""
A non linear static analysis is a collection of parameters for defining the static analysis of the structure.
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
As this analysis is nonlinear the stiffness of the structure is updated at each iteration.
"""
mutable struct NonLinearStaticAnalysis{S<:AbstractStructure,R<:Real,LFV<:Vector{R}} <:
               AbstractStaticAnalysis
    "Structure to be analyzed."
    const s::S
    "Structural state."
    const state::StaticState
    "Load factors vector of the analysis."
    const λᵥ::LFV
    "Current load factor step."
    current_step::Int64
end
"Constructor for a non linear analysis with load factors, optional initial step and initial state."
function NonLinearStaticAnalysis(s::S, λᵥ::LFV;
                                 initial_state::StaticState=StaticState(s),
                                 initial_step::Int=1) where {S<:AbstractStructure,
                                                             LFV<:AbstractVector{<:Real}}
    !(1 ≤ initial_step ≤ length(λᵥ)) &&
        throw(ArgumentError("initial_step must be in [1, $(length(λᵥ))] but is: $initial_step."))
    NonLinearStaticAnalysis(s, initial_state, λᵥ, initial_step)
end

"Constructor for non linear static analysis given a final time (or load factor) and the number of steps."
function NonLinearStaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, initial_step::Int=1)
    t₀ = t₁ / NSTEPS
    λᵥ = collect(LinRange(t₀, t₁, NSTEPS))
    NonLinearStaticAnalysis(s, λᵥ; initial_step)
end

function Base.show(io::IO, sa::NonLinearStaticAnalysis)
    println("NonLinearStaticAnalysis for:")
    println("• Current load factor of $(unwrap(sa.current_step)).")
    show(io, sa.s)
    show(io, sa.state)
end

"Solves an non linear static analysis problem with a given solver."
function _solve!(sa::NonLinearStaticAnalysis, alg::AbstractSolver)
    s = structure(sa)
    # Initialize solution.
    sol = StatesSolution(sa, alg)

    # Load factors iteration.
    while !is_done(sa)

        # Reset assembled magnitudes
        reset!(current_iteration(sa))

        # Computes external forces
        external_forces(current_state(sa)) .= 0
        @debugtime "Assemble external forces" apply!(sa, load_bcs(boundary_conditions(s)))

        # Displacements iteration.
        while isconverged!(current_iteration(sa), tolerances(alg)) isa NotConvergedYet
            # Compute residual forces and tangent matrix.
            @debugtime "Assemble internal forces" assemble!(s, sa)

            # Increment structure displacements `U = U + ΔU`.
            @debugtime "Step" step!(sa, alg)
        end
        # Save current state.
        @debugtime "Save current state" push!(sol, current_state(sa))

        # Increment the time or load factor step.
        @debugtime "Next step" next!(sa)
    end
    sol
end

"Computes ΔU for solving the non linear static analysis with a Newton Raphson method."
function step!(sa::NonLinearStaticAnalysis, ::NewtonRaphson)
    # Extract state info
    state = current_state(sa)
    free_dofs_idx = free_dofs(state)

    # Compute Δu
    r = residual_forces!(state)
    K = tangent_matrix(state)[free_dofs_idx, free_dofs_idx]
    ΔU = Δ_displacements(state)
    cg!(ΔU, K, r)

    # Compute norms
    norm_ΔU = norm(ΔU)
    rel_norm_ΔU = norm_ΔU / norm(displacements(state))
    norm_r = norm(r)
    rel_norm_r = norm_r / norm(external_forces(state))

    # Update displacements into the state.
    state.Uᵏ[free_dofs_idx] .+= ΔU

    # Update iteration
    update!(current_iteration(sa), norm_ΔU, rel_norm_ΔU, norm_r, rel_norm_r)
end

end # module
