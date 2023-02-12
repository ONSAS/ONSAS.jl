"""
Module defining structural solvers that can be used to solved different analyses. 
"""
module StructuralSolvers

export solve
export ConvergenceSettings, step_size, tolerances, step!, AbstractSolution

""" ConvergenceSettings struct.
Facilitates the process of defining and checking numerical convergence. 
### Fields:
- `rel_disp_tol`    -- Relative displacement tolerance.
- `rel_force_tol`   -- Relative residual force tolerance.
- `max_iter`        -- Maximum number of iterations.
"""
Base.@kwdef struct ConvergenceSettings
    rel_disp_tol::Float64 = 1e-6
    rel_force_tol::Float64 = 1e-6
    max_iter::Int = 20
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
function step!(alg::AbstractSolver, analysis::A, args...; kwargs...) where {A} end

include("./../algorithms/NR.jl")

#==========#
# Solutions
#==========#

"""Abstract supertype that holds the solution of a numerical integration."""
abstract type AbstractSolution end

algorithm(sol::AbstractSolution) = sol.alg
displacements(sol::AbstractSolution) = sol.U

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
    initialized_analysis = init(analysis, alg, args...; kwargs...)
    return _solve(initialized_analysis, alg, args...; kwargs...)
end

"Internal solve function to be overloaded by each analysis."
function _solve(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A} end

"Returns the initialized analysis."
function init(analysis::A, alg::AbstractSolver, args...; kwargs...) where {A}
    return analysis
end

end # module