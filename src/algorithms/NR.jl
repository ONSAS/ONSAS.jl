using ..StructuralSolvers: AbstractSolver, ConvergenceSettings

import ..StructuralSolvers: tolerances

export NewtonRaphson

""" Newton-Raphson solver struct.

### Fields:
- `tol` -- Numerical tolerances
"""

struct NewtonRaphson <: AbstractSolver
    tol::ConvergenceSettings
end

# Default tolerances
NewtonRaphson() = NewtonRaphson(ConvergenceSettings())

"Returns NR tolerances settings"
tolerances(alg::NewtonRaphson) = alg.tol