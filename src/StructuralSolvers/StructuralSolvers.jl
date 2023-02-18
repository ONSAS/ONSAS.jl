"""
Module defining structural solvers that can be used to solved different analyses. 
Each solver consists of a data type with a convergence criterion and a the iteration status. 
A _step! method is used to perform a single iteration step.
"""
module StructuralSolvers

using LinearAlgebra: norm

export AbstractConvergenceCriterion, ResidualForceCriterion, ΔUCriterion,
    MaxIterCriterion, ΔU_and_ResidualForce_Criteria, MaxIterCriterion, NotConvergedYet

export ConvergenceSettings, residual_forces_tol, displacement_tol, max_iter_tol
export ResidualsIterationStep, iter, criterion, _reset!, isconverged!, _update!
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

"Returns residual forces tolerance set in the `ConvergenceSettings` `tol`."
residual_forces_tol(tols::ConvergenceSettings) = tols.rel_res_force_tol

"Returns displacements tolerance set in the `ConvergenceSettings` `tol`."
displacement_tol(tols::ConvergenceSettings) = tols.rel_U_tol

"Returns the maximum number of iterations set in the `ConvergenceSettings` `tol`."
max_iter_tol(tols::ConvergenceSettings) = tols.max_iter

""" Abstract supertype for all structural solvers."""

""" Abstract supertype for all convergence criterion."""
abstract type AbstractConvergenceCriterion end

""" `ResidualForceCriterion` convergence criterion. """
struct ResidualForceCriterion <: AbstractConvergenceCriterion end

""" `ΔUCriterion` indicates displacements increment convergence criterion. """
struct ΔUCriterion <: AbstractConvergenceCriterion end

""" `ΔU_and_ResidualForce_Criteria` convergence criterion indicates that both ΔU and residual forces
converged . """
struct ΔU_and_ResidualForce_Criteria <: AbstractConvergenceCriterion end

""" `MaxIterCriterion` criteria indicates that the maximum number of iterations has been reached. """
struct MaxIterCriterion <: AbstractConvergenceCriterion end

""" `NotConvergedYet` indicates that the current iteration has not converged. """
struct NotConvergedYet <: AbstractConvergenceCriterion end


""" ResidualsIterationStep struct.
Stores the convergence information of at the current iteration step.
### Fields:
`Δu_norm`   -- Norm of the displacement increment.
`Δr_norm`   -- Norm of the residual force increment.
`Δu_rel`    -- Relative norm of the displacement increment.
`Δr_rel`    -- Relative norm of the residual force increment.
`n_iter`    -- Current iteration number.
"""
Base.@kwdef mutable struct ResidualsIterationStep{T}
    ΔU_norm::T = 1e12
    Δr_norm::T = 1e12
    ΔU_rel::T = 1e12
    Δr_rel::T = 1e12
    iter::Int = 0
    criterion::AbstractConvergenceCriterion = NotConvergedYet()
end

"Increments a `ResidualsIterationStep` `i_step` by 1."
_step!(i_step::ResidualsIterationStep) = i_step.iter += 1

"Returns the current iteration number of the `ResidualsIterationStep` `i_step`."
iter(i_step::ResidualsIterationStep) = i_step.iter

"Returns the current criterion of the `ResidualsIterationStep` `i_step`."
criterion(ri_step::ResidualsIterationStep) = ri_step.criterion

"Returns the current absolute and relative residual forces norm of the `ResidualsIterationStep` `i_step`."
residual_forces_tol(ri_step::ResidualsIterationStep) = (ri_step.Δr_rel, ri_step.Δr_norm)

"Returns the current absolute and relative residual forces norm of the `ResidualsIterationStep` `i_step`."
displacement_tol(ri_step::ResidualsIterationStep) = (ri_step.ΔU_rel, ri_step.ΔU_norm)

"Updates displacement residuals in the `ResidualsIterationStep` struct `i_step`."
function _update_U!(ri_step::ResidualsIterationStep, U::AbstractVector, ΔU_rel::AbstractVector)
    ri_step.ΔU_norm = norm(ΔU_rel)
    ri_step.ΔU_rel = ri_step.ΔU_norm / norm(U)
end

"Updates forces residuals in the `ResidualsIterationStep` struct `i_step`."
function _update_r!(ri_step::ResidualsIterationStep, fₑₓₜ::AbstractVector{T}, Δr::AbstractVector{T}) where {T<:Real}
    ri_step.Δr_norm = norm(Δr)
    ri_step.Δr_rel = ri_step.Δr_norm / norm(fₑₓₜ)
end

"Updates the `ResidualsIterationStep` `i_step` current convergence criterion."
_update!(ri_step::ResidualsIterationStep, criterion::AbstractConvergenceCriterion) = ri_step.criterion = criterion

"Sets the `ResidualsIterationStep` `i_step` iteration step to 0."
function _reset!(ri_step::ResidualsIterationStep{T}) where {T<:Real}
    ri_step.ΔU_norm = ri_step.Δr_norm = ri_step.ΔU_rel = ri_step.Δr_rel = 1e14 * ones(T)[1]
    ri_step.iter = 0
    ri_step.criterion = NotConvergedYet()
    return ri_step
end

"Updates the `ResidualsIterationStep` `i_step` iteration step with the current values of the
displacement `U` and its increment `ΔU`, the residual increment `Δr`, the external forces vector `fₑₓₜ`
and checks if the iteration has converged given a `ConvergenceSettings` `cs`."
function _update!(ri_step::ResidualsIterationStep,
    ΔU::AbstractVector, U::AbstractVector,
    Δr::AbstractVector, fₑₓₜ::AbstractVector)

    _update_U!(ri_step, U, ΔU)
    _update_r!(ri_step, fₑₓₜ, Δr)
    _step!(ri_step)

    return ri_step
end


"Updates the `ResidualsIterationStep` `i_step` given a `ConvergenceSettings` `cs`."
function isconverged!(ri_step::ResidualsIterationStep, cs::ConvergenceSettings)

    ΔU_relᵏ, ΔU_nromᵏ = displacement_tol(ri_step)
    Δr_relᵏ, Δr_nromᵏ = residual_forces_tol(ri_step)

    @assert ΔU_relᵏ > 0 "Rsidual displacements must be greater than 0."
    @assert Δr_relᵏ > 0 "Residual forces must be greater than 0."

    ΔU_rel_tol = displacement_tol(cs)
    Δr_rel_tol = residual_forces_tol(cs)
    max_iter = max_iter_tol(cs)

    # Check displacements convergence  
    if ΔU_relᵏ ≤ ΔU_rel_tol && ΔU_relᵏ > 0
        _update!(ri_step, ΔUCriterion())
    end

    # Check residual forces convergence 
    if Δr_relᵏ ≤ Δr_rel_tol && Δr_relᵏ > 0
        _update!(ri_step, ResidualForceCriterion())
    end

    if iter(ri_step) > max_iter
        _update!(ri_step, MaxIterCriterion())
        @warn "Maximum number of iterations was reached."
    end

    if Δr_relᵏ ≤ Δr_rel_tol && ΔU_relᵏ ≤ ΔU_rel_tol || Δr_relᵏ < eps() || ΔU_relᵏ < eps()
        _update!(ri_step, ΔU_and_ResidualForce_Criteria())
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

"Returns the step size set to the `AbstractSolver` `solver` ."
step_size(solver::AbstractSolver) = solver.Δt

"Returns the numerical tolerances set to the `AbstractSolver` `solver`."
tolerances(solver::AbstractSolver) = solver.tol

"Computes a step in time on the `analysis` considering the numerical `AbstractSolver` `solver`."
function _step!(solver::AbstractSolver, analysis::A, args...; kwargs...) where {A} end

include("./NewtonRaphson.jl")

# ===============
# Solve function
# ===============

"""
Solve an structural analysis problem.
### Input
- `analysis` -- structural analysis problem
- `alg`     -- structural algorithm to solve the problem
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

include("./Assembler.jl")

end # module