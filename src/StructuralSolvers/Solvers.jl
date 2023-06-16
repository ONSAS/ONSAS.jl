module Solvers

using Reexport

using ..StructuralSolvers

@reexport import ..StructuralSolvers: tolerances

export NewtonRaphson

"""
Newton-Raphson solver struct.
"""
struct NewtonRaphson <: AbstractSolver
    "Numerical tolerances."
    tol::ConvergenceSettings
end

"Constructor for `NewtonRaphson` solver with default tolerances."
NewtonRaphson() = NewtonRaphson(ConvergenceSettings())

"Return `NewtonRaphson` solver tolerances."
tolerances(sol::NewtonRaphson) = sol.tol

end # module
