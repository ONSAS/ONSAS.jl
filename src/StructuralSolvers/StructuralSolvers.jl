"""
Module defining structural solvers that can be used to solved different analyses. 
Each solver consists of a data type with a convergence criterion and a the iteration status. 
A _step! method is used to perform a single iteration step.
"""
module StructuralSolvers

using LinearAlgebra: norm

using ..Elements

export AbstractConvergenceCriterion, ResidualForceCriterion, ΔUCriterion,
       MaxIterCriterion, ΔU_and_ResidualForce_Criteria, MaxIterCriterion, NotConvergedYet
export ConvergenceSettings, residual_forces_tol, displacement_tol, max_iter_tol
export ResidualsIterationStep, iter, criterion, _reset!, isconverged!, _update!
export AbstractSolver, step_size, tolerances, _step!, solve!, _solve!, solve, reset!
export AbstractSolution

"""
Facilitates the process of defining and checking numerical convergence. 
"""
Base.@kwdef struct ConvergenceSettings
    "Relative displacement tolerance."
    rel_U_tol::Float64 = 1e-6
    "Relative residual force tolerance."
    rel_res_force_tol::Float64 = 1e-6
    "Maximum number of iterations."
    max_iter::Int = 20
end

"Return residual forces tolerance set in the `ConvergenceSettings` `tol`."
residual_forces_tol(tols::ConvergenceSettings) = tols.rel_res_force_tol

"Return displacements tolerance set in the `ConvergenceSettings` `tol`."
displacement_tol(tols::ConvergenceSettings) = tols.rel_U_tol

"Return the maximum number of iterations set in the `ConvergenceSettings` `tol`."
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

"""
Stores the convergence information of at the current iteration step.
"""
Base.@kwdef mutable struct ResidualsIterationStep{T}
    "Norm of the displacement increment."
    ΔU_norm::T = 1e12
    "Norm of the residual force increment."
    Δr_norm::T = 1e12
    "Relative norm of the displacement increment."
    ΔU_rel::T = 1e12
    "Relative norm of the residual force increment."
    Δr_rel::T = 1e12
    "Current iteration number."
    iter::Int = 0
    criterion::AbstractConvergenceCriterion = NotConvergedYet()
end

"Increments a `ResidualsIterationStep` `i_step` by 1."
_step!(i_step::ResidualsIterationStep) = i_step.iter += 1

"Return the current iteration number of the `ResidualsIterationStep` `i_step`."
iter(i_step::ResidualsIterationStep) = i_step.iter

"Return the current criterion of the `ResidualsIterationStep` `i_step`."
criterion(ri_step::ResidualsIterationStep) = ri_step.criterion

"Return the current absolute and relative residual forces norm of the `ResidualsIterationStep` `i_step`."
residual_forces_tol(ri_step::ResidualsIterationStep) = (ri_step.Δr_rel, ri_step.Δr_norm)

"Return the current absolute and relative residual forces norm of the `ResidualsIterationStep` `i_step`."
displacement_tol(ri_step::ResidualsIterationStep) = (ri_step.ΔU_rel, ri_step.ΔU_norm)

"Sets the `ResidualsIterationStep` `i_step` iteration step to 0."
function _reset!(ri_step::ResidualsIterationStep{T}) where {T<:Real}
    ri_step.ΔU_norm = ri_step.Δr_norm = ri_step.ΔU_rel = ri_step.Δr_rel = 1e14 * ones(T)[1]
    ri_step.iter = 0
    ri_step.criterion = NotConvergedYet()
    return ri_step
end

"Sets the `ResidualsIterationStep` `i_step` iteration step with nothing."
function _reset!(ri_step::ResidualsIterationStep{<:Nothing})
    ri_step.iter = 0
    ri_step.criterion = ResidualForceCriterion()
    return ri_step
end

"Updates the `ResidualsIterationStep` `i_step` current convergence criterion."
function _update!(ri_step::ResidualsIterationStep, criterion::AbstractConvergenceCriterion)
    return ri_step.criterion = criterion
end

"Updates the `ResidualsIterationStep` `i_step` iteration step with the current values of the displacement residuals
`ΔU_norm` and `ΔU_rel`, the norm of the residual forces increment `Δr_norm`, and the relative `Δr_rel`."
function _update!(ri_step::ResidualsIterationStep, ΔU_norm::Real, ΔU_rel::Real, Δr_norm::Real,
                  Δr_rel::Real)
    ri_step.ΔU_norm = ΔU_norm
    ri_step.ΔU_rel = ΔU_rel
    ri_step.Δr_norm = Δr_norm
    ri_step.Δr_rel = Δr_rel

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

    if Δr_relᵏ ≤ Δr_rel_tol && ΔU_relᵏ ≤ ΔU_rel_tol || Δr_nromᵏ < eps() || ΔU_nromᵏ < eps()
        _update!(ri_step, ΔU_and_ResidualForce_Criteria())
        return ΔU_and_ResidualForce_Criteria()
    end

    return NotConvergedYet()
end

#==========#
# Solvers
#==========#

"""
Abstract supertype for all direct integration methods.
"""
abstract type AbstractSolver end

"Return the step size set to the `AbstractSolver` `solver` ."
step_size(solver::AbstractSolver) = solver.Δt

"Return the numerical tolerances set to the `AbstractSolver` `solver`."
tolerances(solver::AbstractSolver) = solver.tol

"Computes a step in time on the `analysis` considering the numerical `AbstractSolver` `solver`."
function _step!(solver::AbstractSolver, analysis::A) where {A} end

# ===============
# Solve function
# ===============

"""
Solve a structural analysis problem with the given solver.
"""
function solve(analysis::A, alg::AbstractSolver=nothing, args...; kwargs...) where {A}
    # FIXME Errors copying the mesh struct.
    analysis = deepcopy(analysis)
    reset!(analysis)
    solve!(analysis, alg, args...; kwargs...)
end

"""
Solve a structural analysis problem with the given solver.
This function mutates the state defined in the `analysis` problem; use [`solve`](@ref) to avoid mutation.
For linear analysis problems, `alg` doesn't need to be provided.

Returns a solution structure holding the result and the algorithm used to obtain it.
"""
function solve!(analysis::A, alg::Union{AbstractSolver,Nothing}=nothing, args...;
                kwargs...) where {A}
    # TODO Dispatch on `nothing` only for linear problems.
    if isnothing(alg)
        _solve!(analysis, args...; kwargs...)
    else
        initialized_analysis = _init(analysis, alg, args...; kwargs...)
        _solve!(initialized_analysis, alg, args...; kwargs...)
    end
end

"Internal solve function to be overloaded by each analysis."
function _solve!(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A} end

"Return the initialized `analysis`. By default, it Return the same analysis."
_init(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A} = analysis

"Resets the analysis to the state before starting a new assembly."
function reset! end

end # module
