"""
Module defining structural solvers that can be used to solved different analyses. 
"""
module StructuralSolvers

using Reexport: @reexport
@reexport import ..StructuralModel: displacements
@reexport import ..StructuralAnalyses: AbstractStructuralAnalysis
@reexport import ..Utils: solve

export ConvergenceSettings, step_size, tolerances, NewtonRaphson, AbstractSolution, StaticSolution


const DEFAULT_RELATIVE_DISP_TOLERANCE = 1e-6
const DEFAULT_RELATIVE_FORCE_TOLERANCE = 1e-6
const DEFAULT_MAX_ITER = 20

""" ConvergenceSettings struct.
Facilitates the process of defining and checking numerical convergence. 
### Fields:
- `rel_disp_tol`    -- Relative displacement tolerance.
- `rel_force_tol`   -- Relative residual force tolerance.
- `max_iter`        -- Maximum number of iterations.
"""
Base.@kwdef struct ConvergenceSettings
    rel_disp_tol::Number = DEFAULT_RELATIVE_DISP_TOLERANCE
    rel_force_tol::Number = DEFAULT_RELATIVE_FORCE_TOLERANCE
    max_iter::Integer = DEFAULT_MAX_ITER
end

#==========#
# Solvers
#==========#

"""
    AbstractSolver
Abstract supertype for all direct integration methods.
"""
abstract type AbstractSolver end
"Returns the step size of the solver."
step_size(alg::AbstractSolver) = alg.Δt
"Returns the numerical tolerance set to the solver."
tolerances(alg::AbstractSolver) = alg.tol

"Computes the solver step for a given analysis."
function step!(alg::AbstractSolver, a::AbstractStructuralAnalysis, args...; kwargs...) end

include("./../algorithms/NR.jl")

#==========#
# Solutions
#==========#

"""
    AbstractSolution
Abstract supertype that holds the solution of a numerical integration.
"""
abstract type AbstractSolution end

algorithm(sol::AbstractSolution) = sol.alg
displacements(sol::AbstractSolution) = sol.U

"""
Static solution struct.
### Fields
- `alg` -- Algorithm used in the integration
- `U`   -- Displacements
- `t`   -- Vector of time values
- `vars`-- Other relevant variables are stored in a dictionary
"""
struct StaticSolution{T<:AbstractSolver,UT,ST} <: AbstractSolution
    alg::T
    U::UT
    t::ST
    vars::Dict
end

# constructor with missing fields
StaticSolution(alg, U, t) = StaticSolution(alg, U, t, Dict())

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
function solve(a::AbstractStructuralAnalysis, alg::AbstractSolver, args...; kwargs...)
    initialized_analysis = init(a, alg, args...; kwargs...)
    step!(alg, initialized_analysis, args...; kwargs...)
    return _solve(initialized_analysis, alg, args...; kwargs...)
end

"Returns the initialized analysis."
function init(a::AbstractStructuralAnalysis, alg::AbstractSolver, args...; kwargs...)
    return a
end

end # module