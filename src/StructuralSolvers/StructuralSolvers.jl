"""
Module defining structural solvers that can be used to solved different analyses.
Each solver consists of a data type with a convergence criterion and a the iteration status.
A step! method is used to perform a single iteration step.
"""
module StructuralSolvers

using LinearAlgebra: norm
using Reexport
using LinearSolve
using IterativeSolvers

using ..StructuralAnalyses
using ..Entities

@reexport import ..Assemblers: reset!
@reexport import ..StructuralAnalyses: tangent_matrix
@reexport import ..CommonSolve: solve!, init

export AbstractConvergenceCriterion, ResidualForceCriterion, ΔUCriterion,
       MaxIterCriterion, ΔU_and_ResidualForce_Criteria, MaxIterCriterion, NotConvergedYet,
       ConvergenceSettings, residual_forces_tol, displacement_tol, max_iter_tol,
       ResidualsIterationStep, iter, criterion, isconverged!,
       AbstractSolver, step_size, tolerances, step!, solve, solve!, _solve!,
       AbstractSolution, iterations, update!, next!

const INITIAL_Δ = 1e12
const DEFAULT_LINEAR_SOLVER = IterativeSolversJL_CG
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

"Show convergence settings."
function Base.show(io::IO, cs::ConvergenceSettings)
    println("• ||ΔUᵏ||/||Uᵏ||| ≤ : $(residual_forces_tol(cs))")
    println("• ||ΔRᵏ||/||Fₑₓₜ||| $(displacement_tol(cs))")
    println("• iter k ≤: $(max_iter_tol(cs))")
end

"Return residual forces tolerance set."
residual_forces_tol(tols::ConvergenceSettings) = tols.rel_res_force_tol

"Return displacements tolerance set."
displacement_tol(tols::ConvergenceSettings) = tols.rel_U_tol

"Return the maximum number of iterations set."
max_iter_tol(tols::ConvergenceSettings) = tols.max_iter

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
    ΔU_norm::T = INITIAL_Δ
    "Norm of the residual force increment."
    Δr_norm::T = INITIAL_Δ
    "Relative norm of the displacement increment."
    ΔU_rel::T = INITIAL_Δ
    "Relative norm of the residual force increment."
    Δr_rel::T = INITIAL_Δ
    "Current iteration number."
    iter::Int = 0
    criterion::AbstractConvergenceCriterion = NotConvergedYet()
end

"Increments a `ResidualsIterationStep` `i_step` by 1."
step!(i_step::ResidualsIterationStep) = i_step.iter += 1

"Return the iterations done so far."
iterations(i_step::ResidualsIterationStep) = i_step.iter

"Return the current convergence criterion."
criterion(ri_step::ResidualsIterationStep) = ri_step.criterion

"Return the current absolute and relative residual forces norm."
residual_forces_tol(ri_step::ResidualsIterationStep) = (ri_step.Δr_rel, ri_step.Δr_norm)

"Return the current absolute and relative residual forces norm."
displacement_tol(ri_step::ResidualsIterationStep) = (ri_step.ΔU_rel, ri_step.ΔU_norm)

"Sets the iteration step to 0."
function reset!(ri_step::ResidualsIterationStep{T}) where {T<:Real}
    ri_step.ΔU_norm = ri_step.Δr_norm = ri_step.ΔU_rel = ri_step.Δr_rel = INITIAL_Δ * ones(T)[1]
    ri_step.iter = 0
    ri_step.criterion = NotConvergedYet()
    ri_step
end

"Sets the iteration step with nothing."
function reset!(ri_step::ResidualsIterationStep{<:Nothing})
    ri_step.iter = 0
    ri_step.criterion = ResidualForceCriterion()
    ri_step
end

"Updates the current convergence criterion."
function update!(ri_step::ResidualsIterationStep, criterion::AbstractConvergenceCriterion)
    ri_step.criterion = criterion
end

"Updates the iteration step with the current values of the displacement and forces residuals."
function update!(ri_step::ResidualsIterationStep, ΔU_norm::Real, ΔU_rel::Real, Δr_norm::Real,
                 Δr_rel::Real)
    ri_step.ΔU_norm = ΔU_norm
    ri_step.ΔU_rel = ΔU_rel
    ri_step.Δr_norm = Δr_norm
    ri_step.Δr_rel = Δr_rel
    step!(ri_step)
    ri_step
end

"Updates the convergence criteria."
function isconverged!(ri_step::ResidualsIterationStep, cs::ConvergenceSettings)
    ΔU_relᵏ, ΔU_normᵏ = displacement_tol(ri_step)
    Δr_relᵏ, Δr_normᵏ = residual_forces_tol(ri_step)

    @assert ΔU_relᵏ > 0 "Residual displacements norm must be greater than 0."
    @assert Δr_relᵏ > 0 "Residual forces norm must be greater than 0."

    ΔU_rel_tol = displacement_tol(cs)
    Δr_rel_tol = residual_forces_tol(cs)
    max_iter = max_iter_tol(cs)

    criterion = if Δr_relᵏ ≤ Δr_rel_tol
        ResidualForceCriterion()
        # Check residual forces convergence
    elseif ΔU_relᵏ ≤ ΔU_rel_tol
        ΔUCriterion()
    elseif iterations(ri_step) > max_iter
        @warn "Maximum number of iterations was reached."
        MaxIterCriterion()
    else
        NotConvergedYet()
    end
    update!(ri_step, criterion)
end

#==========#
# Solvers
#==========#

"""
Abstract supertype for all direct integration methods.
"""
abstract type AbstractSolver end

"Return the step size."
step_size(solver::AbstractSolver) = solver.Δt

"Return the numerical tolerances."
tolerances(solver::AbstractSolver) = solver.tol

"Computes a step in time on the `analysis` considering the numerical `AbstractSolver` `solver`."
function step!(solver::AbstractSolver, analysis::A) where {A} end

"Increment the time step given of a structural analysis. Dispatch is done for different
solvers."
next!(a::AbstractStructuralAnalysis, solver::AbstractSolver) = a.t += time_step(a) # TODO Define `time_step` fallback.

"Return system tangent matrix in the structural state given a solver."
function tangent_matrix(st::AbstractStructuralState, alg::AbstractSolver) end

# ===============
# Solve function
# ===============

"""
Solve a structural analysis problem with the given solver, returning a solution structure
which holds the result and the algorithm used to obtain it.

This function mutates the state defined in the analysis problem.
Use [`solve`](@ref) to avoid mutation.
For linear analysis problems, the algorithm doesn't need to be provided.
Also a linar solver form the LinearSolve.jl package can by provided. By default the wrapper to IterativeSolvers Conjugatge Gradient is being used.
"""
function solve!(problem::AbstractStructuralAnalysis, solver::Union{AbstractSolver,Nothing}=nothing,
                linear_solve::SciMLBase.AbstractLinearAlgorithm=DEFAULT_LINEAR_SOLVER())
    _solve!(problem, solver, linear_solve)
end
solve!(pair::Tuple{AbstractStructuralAnalysis,AbstractSolver}) = solve!(pair.first, pair.last)

"Copy the structure and optionally reset the state."
function init(a::AbstractStructuralAnalysis, solver::AbstractSolver; reset::Bool=true)
    acopy = deepcopy(a)
    (reset ? reset!(acopy) : acopy, solver)
end

"Return default linear solver tolerances"
function _default_linear_solver_tolerances(A::AbstractMatrix{<:Real}, b::Vector{<:Real})
    abstol = zero(real(eltype(b)))
    reltol = sqrt(eps(real(eltype(b))))
    maxiter = length(b)
    abstol, reltol, maxiter
end

"Internal solve function to be overloaded by each analysis."
function _solve!(problem::AbstractStructuralAnalysis, solver::AbstractSolver, args...; kwargs...) end

end # module
