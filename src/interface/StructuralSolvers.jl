"""
Module defining structural solvers that can be used to solved different analyses. 
"""
module StructuralSolvers

using LinearAlgebra: norm

export AbstractConvergenceCriterion, ResidualForceCriterion, ΔUCriterion,
    MaxIterCriterion, ΔU_and_ResidualForce_Criteria, MaxIterCriterion, NotConvergedYet, ConvergenceSettings, _has_converged!
export IterationStep, _reset!, _has_converged!, _update!
export AbstractSolver, step_size, tolerances, _step!, solve, _solve

""" ConvergenceSettings struct.
Facilitates the process of defining and checking numerical convergence. 
### Fields:
- `rel_disp_tol`    -- Relative displacement tolerance.
- `rel_force_tol`   -- Relative residual force tolerance.
- `max_iter`        -- Maximum number of iterations.
"""
Base.@kwdef struct ConvergenceSettings
    rel_U_tol::Float64 = 1e-6
    rel_res_force_tol::Float64 = 1e-6
    max_iter::Int = 20
end

""" Abstract supertype for all convergence criteria."""
abstract type AbstractConvergenceCriterion end

""" Residual force convergence criterion. """
struct ResidualForceCriterion <: AbstractConvergenceCriterion end

""" Displacements increment convergence criterion. """
struct ΔUCriterion <: AbstractConvergenceCriterion end

""" Both criterion convergence. """
struct ΔU_and_ResidualForce_Criteria <: AbstractConvergenceCriterion end

""" Maximum number of iterations convergence criterion. """
struct MaxIterCriterion <: AbstractConvergenceCriterion end

""" Not converged criterion. """
struct NotConvergedYet <: AbstractConvergenceCriterion end


""" IterationStep struct.
Stores the information of a single iteration step.
### Fields:
TODO
"""
Base.@kwdef mutable struct IterationStep{T}
    Δu_norm::T = 0.0
    Δr_norm::T = 0.0
    Δu_rel::T = 0.0
    Δr_rel::T = 0.0
    iter::Int = 0
    criterion::AbstractConvergenceCriterion = NotConvergedYet()
end

"Rests the iteration step."
function _reset!(i_step::IterationStep{T}) where {T<:Real}
    i_step.Δu_norm = i_step.Δr_norm = i_step.Δu_rel = i_step.Δr_rel = zeros(T)[1]
    i_step.iter = 0
    i_step.criterion = NotConvergedYet()
    return nothing
end

"Updates the iteration step."
function _update!(i_step::IterationStep,
    ΔU::AbstractVector, U::AbstractVector,
    r::AbstractVector, fₑₓₜ::AbstractVector,
    cs::ConvergenceSettings)

    i_step.Δu_norm = norm(ΔU)
    i_step.Δu_rel = i_step.Δu_norm / norm(U)

    i_step.Δr_norm = norm(r)
    i_step.Δr_rel = i_step.Δr_norm / norm(fₑₓₜ)

    i_step.iter += 1

    _has_converged!(i_step, cs)

    return i_step
end


"Checks if a given iteration has converged."
function _has_converged!(i_step::IterationStep, cs::ConvergenceSettings)

    Δu_rel = i_step.Δu_rel
    Δr_rel = i_step.Δr_rel

    if Δu_rel ≤ cs.rel_U_tol && Δu_rel > 0
        i_step.criterion = ΔUCriterion()
    end

    if Δr_rel ≤ cs.rel_res_force_tol && Δr_rel > 0
        i_step.criterion = ResidualForceCriterion()
    end

    if i_step.iter > cs.max_iter
        i_step.criterion = MaxIterCriterion()
        @warn "Maximum number of iterations reached."
    end

    if Δu_rel ≤ cs.rel_U_tol && Δr_rel ≤ cs.rel_res_force_tol && Δu_rel > 0 && Δr_rel > 0
        i_step.criterion = ΔU_and_ResidualForce_Criteria()
        return true
    end

    return false
end

#==========#
# Solvers
#==========#

"""
Abstract supertype for all direct integration methods.
"""
abstract type AbstractSolver end

"Returns the step size of the solver."
step_size(alg::AbstractSolver) = alg.Δt

"Returns the numerical tolerance set to the solver."
tolerances(alg::AbstractSolver) = alg.tol

"Computes solver step for a given analysis."
function _step!(alg::AbstractSolver, analysis::A, args...; kwargs...) where {A} end

include("./../algorithms/NR.jl")

# ===============
# Solve function
# ===============

"""
Solve an structural analysis problem.
### Input
- `a` -- structural analysis problem
- `alg` -- structural algorithm to solve the problem
### Output
A solution structure (`AbstractSolution`) that holds the result and the algorithm used
to obtain it.
"""
function solve(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A}
    initialized_analysis = _init(analysis, alg, args...; kwargs...)
    return _solve(initialized_analysis, alg, args...; kwargs...)
end

"Internal solve function to be overloaded by each analysis."
function _solve(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A} end

"Returns the initialized analysis."
function _init(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A}
    return analysis
end

end # module