module NonLinearStaticAnalyses

using LinearAlgebra: norm
using IterativeSolvers: cg

using ....Utils: ScalarWrapper, @debugtime
using ..StaticAnalyses
using ...StructuralSolvers: AbstractSolver, NewtonRaphson, StatesSolution, tolerances
using ...StructuralModel: AbstractStructure

import ...StructuralSolvers: _solve!, _step!

export NonLinearStaticAnalysis

""" NonLinearStaticAnalysis struct.
A `NonLinearStaticAnalysis` is a collection of parameters for defining the static analysis of the structure. 
In the static analysis, the structure is analyzed at a given load factor (this variable is analog to time).
As this analysis is nonlinear the stiffness of the structure is updated at each iteration. 
### Fields:
- `s`             -- stores the structure to be analyzed.
- `state`         -- stores the structural state.
- `λᵥ`            -- stores the load factors vector of the analysis
- `current_step`  -- stores the current load factor step
"""
struct NonLinearStaticAnalysis{S<:AbstractStructure,LFV<:AbstractVector{<:Real}} <: AbstractStaticAnalysis
    s::S
    state::StaticState
    λᵥ::LFV
    current_step::ScalarWrapper{Int}
    function NonLinearStaticAnalysis(s::S, λᵥ::LFV; initial_step::Int=1) where {S<:AbstractStructure,LFV<:AbstractVector{<:Real}}
        new{S,LFV}(s, StaticState(s), λᵥ, ScalarWrapper(initial_step))
    end
end

"Constructor for `NonLinearStaticAnalysis` given a final time (or load factor) `t₁` and the number of steps `NSTEPS`."
function NonLinearStaticAnalysis(s::AbstractStructure, t₁::Real=1.0; NSTEPS=10, initial_step::Int=1)
    t₀ = t₁ / NSTEPS
    λᵥ = LinRange(t₀, t₁, NSTEPS) |> collect
    NonLinearStaticAnalysis(s, λᵥ, initial_step=initial_step)
end

"Solves an `NonLinearStaticAnalysis` `sa` with an AbstractSolver `alg`."
function _solve!(sa::NonLinearStaticAnalysis, alg::AbstractSolver)
    s = structure(sa)
    # Initialize solution.
    sol = StatesSolution(sa, alg)

    # Load factors iteration.
    while !is_done(sa)

        # Sets Δu, ΔR and relatives norms to zero 
        _reset!(current_iteration(sa))

        # Computes external forces
        _apply!(sa, load_bcs(s))

        # Displacements iteration.
        while isconverged!(current_iteration(sa), tolerances(alg)) isa NotConvergedYet
            # Compute residual forces and tangent matrix.
            @debugtime "Assemble" _assemble!(s, sa)

            # Increment structure displacements `U = U + ΔU`.
            @debugtime "Step" _step!(sa, alg)
        end
        # Save current state.
        @debugtime "Save current state" push!(sol, current_state(sa))

        # Increment the time or load factor step.
        @debugtime "Next step" _next!(sa)
    end
    return sol
end

"Computes ΔU for solving the `NonLinearStaticAnalysis` `sa` with a `NewtonRaphson` method."
function _step!(sa::NonLinearStaticAnalysis, ::NewtonRaphson)

    # Extract state info
    c_state = current_state(sa)
    f_dofs_indexes = index.(free_dofs(c_state))

    # Compute Δu
    r = residual_forces(c_state)
    K = view(tangent_matrix(c_state), f_dofs_indexes, f_dofs_indexes)
    ΔU = cg(K, r)

    # Compute norms
    norm_ΔU = norm(ΔU)
    rel_norm_ΔU = norm_ΔU / norm(displacements(c_state))
    norm_r = norm(r)
    rel_norm_r = norm_r / norm(external_forces(c_state))

    # Update displacements into the state
    _update!(c_state, ΔU)

    # Update iteration 
    _update!(current_iteration(sa), norm_ΔU, rel_norm_ΔU, norm_r, rel_norm_r)

end

end #end module
