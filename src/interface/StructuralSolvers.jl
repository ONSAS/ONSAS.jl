"""
Module defining structural solvers that can be used to solved different analyses. 
"""
module StructuralAnalyses

using Reexport: @reexport
@reexport import ..StructuralModel: displacements
@reexport import ..StructuralAnalyses: AbstractStructuralAnalysis

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
step_size(alg::AbstractSolver) = alg.Δt
tolerances(alg::AbstractSolver) = alg.tol


""" Newton-Raphson solver struct.

### Fields:
- `Δt`  -- Step size
- `tol` -- Numerical tolerances
"""

struct NewtonRaphson <: AbstractSolver
    Δt::Number
    tol::ConvergenceSettings
end

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
    Solution{T<:AbstractSolver, UT, VT, AT} <: AbstractSolution
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
Solution(alg, U, t) = Solution(alg, U, nothing, nothing, t)

""" Returns the ambient dimension of the state space of the solution."""
dim(sol::Solution) = length(first(sol.U))

"""Returns the vector of displacements of the given solution along coordinate `i`."""
function displacements(sol::Solution, i::Int)
    1 ≤ i ≤ dim(sol) || throw(ArgumentError("expected the coordinate to be between 1 and $(dim(sol)), got $i"))
    U = displacements(sol)
    return [u[i] for u in U]
end

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
function solve(a::AbstractStructuralAnalyses, alg::AbstractSolver, args...; kwargs...)
    initialized_analysis = init(ivp, alg, args...; kwargs...)
    return _solve(initialized_analysis, alg, args...; kwargs...)
end

end # module