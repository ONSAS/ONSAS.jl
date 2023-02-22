using ..StructuralSolvers: AbstractSolver, ConvergenceSettings

import ..StructuralSolvers: tolerances

export NewtonRaphson

""" Newton-Raphson solver struct.
### Fields:
- `tols` -- Numerical tolerances
"""
struct NewtonRaphson <: AbstractSolver
    tol::ConvergenceSettings
end

"Constructor for `NewtonRaphson` solver with default tolerances."
NewtonRaphson() = NewtonRaphson(ConvergenceSettings())

"Returns `NewtonRaphson` solver tolerances."
tolerances(sol::NewtonRaphson) = sol.tol